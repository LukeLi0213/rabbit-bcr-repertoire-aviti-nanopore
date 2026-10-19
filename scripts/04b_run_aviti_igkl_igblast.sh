#!/bin/bash
# ============================================================
# 04b_run_aviti_igkl_igblast.sh
# 接續04a_run_aviti_igkl_demux.sh：對每個chunk fasta跑igblastn。
# 架構完全比照02d（chunk+子批次+checkpoint），但「驗證+重跑」機制改用
# round-based capped retry（igblast_round_retry.py，深度上限5層，超過
# 上限的reads收集起來，等這一輪全部batch跑完後，下一輪重新從頭批次
# 處理），取代02d原本的bash遞迴切半版本——10K測試證實IGL在遞迴式
# 重跑下會深切到10層以上、疑似觸發background process異常終止；
# round-based capped retry已在IGK/IGL各10K reads（r1only+r2only）
# 驗證可穩定收斂（IGK 0重試；IGL最多5輪收斂到個位數不可處理reads）。
#
# D-gene佔位符：IGK/IGL沒有D基因，但igblastn 1.22.0沒給-germline_db_D
# 會報錯，借用該chain自己的J基因fasta當佔位符（已驗證：d_call正確
# 維持空白，不產生假陽性）。
#
# 用法：bash 04b_run_aviti_igkl_igblast.sh <rabbit> <IGK|IGL> <r1only|r2only>
# 中斷後重跑：已完成的chunk/子批次會被跳過（.done標記）
# ============================================================
set -e
RABBIT="$1"; CHAIN="$2"; SIDE="$3"
[ -z "$RABBIT" ] || [ -z "$CHAIN" ] || [ -z "$SIDE" ] && {
    echo "用法：bash 04b_run_aviti_igkl_igblast.sh <rabbit> <IGK|IGL> <r1only|r2only>"; exit 1; }

BASE_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
PYTHON="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python"
ROUND_RETRY="${BASE_DIR}/script/igblast_round_retry.py"
IGBLAST_BATCH=25000

if [ "$SIDE" = "r1only" ]; then SIDE_TAG="R1only"; elif [ "$SIDE" = "r2only" ]; then SIDE_TAG="R2only"; else
    echo "SIDE必須是r1only或r2only"; exit 1; fi

OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_${SIDE_TAG}"
CHUNK_DIR="${OUT_DIR}/chunks"
SAMPLE="15651-LT-${RABBIT}_AVITI_${CHAIN}_${SIDE_TAG}"
FINAL_AIRR="${OUT_DIR}/${SAMPLE}.airr.tsv"
STATS="${OUT_DIR}/${SAMPLE}_qc_stats.json"

GV="${IMGT}/${CHAIN}V.fasta"
GJ="${IMGT}/${CHAIN}J.fasta"
GD="${IMGT}/${CHAIN}J.fasta"   # IGK/IGL無D基因，借用J基因檔案佔位

GREEN='\033[0;32m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }

split_fasta() {
    awk -v batch="$3" -v prefix="$2" '
    /^>/ { if (n % batch == 0) { if (fi > 0) close(outfile); outfile = sprintf("%s_%04d.fasta", prefix, fi); fi++ } n++ }
    { print > outfile }' "$1"
}

# 2026-09-30新增：最終合併檔已存在就直接跳過，不重新執行（避免chunk中間
# 檔案被清理後，重跑時「合併所有chunk」步驟報錯或覆蓋掉已完成的結果）。
if [ -s "$FINAL_AIRR" ]; then
    info "最終合併檔已存在（${FINAL_AIRR}），視為此rabbit/chain/side已完成，略過重新執行。"
    info "若確定要強制重跑，請先手動搬移或刪除該檔案，並確認chunks/底下的中間檔案仍完整存在。"
    exit 0
fi

[ -f "$ROUND_RETRY" ] || error "找不到 ${ROUND_RETRY}"
[ -d "$CHUNK_DIR" ] || error "找不到 ${CHUNK_DIR}，請先執行 04a_run_aviti_igkl_demux.sh ${RABBIT} ${CHAIN}"
shopt -s nullglob
FASTAS=("${CHUNK_DIR}"/chunk_*.fasta)
[ ${#FASTAS[@]} -gt 0 ] || error "找不到chunk fasta"
info "Rabbit ${RABBIT} ${CHAIN} ${SIDE}：找到 ${#FASTAS[@]} 個 chunk fasta"

for FASTA in "${FASTAS[@]}"; do
    IDX=$(basename "$FASTA" .fasta | sed 's/chunk_//')
    AIRR="${CHUNK_DIR}/chunk_${IDX}.airr.tsv"
    DONE_MARK="${CHUNK_DIR}/chunk_${IDX}.igblast.done"
    SUB_DIR="${CHUNK_DIR}/chunk_${IDX}_sub"

    if [ -f "$DONE_MARK" ]; then
        info "chunk ${IDX} igblast 已完成，跳過"
        continue
    fi

    N_SEQ=$(grep -c "^>" "$FASTA")
    step "IgBLAST chunk ${IDX}（${N_SEQ} 條，切成 ${IGBLAST_BATCH}/批）"
    mkdir -p "$SUB_DIR"
    SUB_FASTAS=("${SUB_DIR}"/sub_*.fasta)
    if [ ${#SUB_FASTAS[@]} -eq 0 ]; then
        split_fasta "$FASTA" "${SUB_DIR}/sub" "$IGBLAST_BATCH"
        SUB_FASTAS=("${SUB_DIR}"/sub_*.fasta)
    fi
    info "  分成 ${#SUB_FASTAS[@]} 個子批次"

    for SUB_FASTA in "${SUB_FASTAS[@]}"; do
        SUB_ID=$(basename "$SUB_FASTA" .fasta)
        SUB_AIRR="${SUB_DIR}/${SUB_ID}.airr.tsv"
        SUB_DONE="${SUB_DIR}/${SUB_ID}.done"
        [ -f "$SUB_DONE" ] && { info "    ${SUB_ID} 已完成，跳過"; continue; }
        SUB_N=$(grep -c "^>" "$SUB_FASTA")
        info "    跑 ${SUB_ID}（${SUB_N} 條，round-based retry）..."
        "$PYTHON" "$ROUND_RETRY" "$SUB_FASTA" "$SUB_AIRR" "$GV" "$GD" "$GJ"
        touch "$SUB_DONE"
    done

    SUB_AIRRS=("${SUB_DIR}"/sub_*.airr.tsv)
    head -1 "${SUB_AIRRS[0]}" > "$AIRR"
    for SUB_AIRR in "${SUB_AIRRS[@]}"; do tail -n +2 "$SUB_AIRR" >> "$AIRR"; done
    rm -rf "$SUB_DIR"
    touch "$DONE_MARK"
    info "  ✓ chunk ${IDX} 完成：$(tail -n +2 "$AIRR" | wc -l | tr -d ' ') rows"
done

step "合併所有 chunk"
AIRRS=("${CHUNK_DIR}"/chunk_*.airr.tsv)
head -1 "${AIRRS[0]}" > "$FINAL_AIRR"
for AIRR in "${AIRRS[@]}"; do tail -n +2 "$AIRR" >> "$FINAL_AIRR"; done
TOTAL_N=$(tail -n +2 "$FINAL_AIRR" | wc -l | tr -d ' ')
info "合併完成：${FINAL_AIRR}（${TOTAL_N} rows）"

"$PYTHON" -c "
import json, csv
n = 0; prod = 0
with open('${FINAL_AIRR}') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        n += 1
        if row.get('productive') == 'T': prod += 1
with open('${STATS}', 'w') as f:
    json.dump({'sample': '${SAMPLE}', 'platform': 'AVITI', 'chain': '${CHAIN}',
               'side': '${SIDE}', 'method': 'round-based capped retry igblast',
               'reads_igblast': n, 'productive_count': prod}, f, indent=2)
print(f'productive={prod/n*100:.1f}%' if n else 'n=0')
"
info "✓ Rabbit ${RABBIT} ${CHAIN} ${SIDE} igblast 全部完成：${FINAL_AIRR}"

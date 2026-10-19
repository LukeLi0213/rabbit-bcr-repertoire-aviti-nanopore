#!/bin/bash
# ============================================================
# 04d_stitch_igblast_igkl.sh
# 接續04c_stitch_qc_igkl.py：對拼接+QC通過的chunk fasta跑第3次igblast，
# 取得完整V(D)J判讀。架構跟04b完全相同（chunk+子批次+checkpoint+
# round-based capped retry），只是輸入/輸出目錄換成{CHAIN}_Stitched。
#
# 用法：bash 04d_stitch_igblast_igkl.sh <rabbit> <IGK|IGL>
# 中斷後重跑：已完成的chunk/子批次會被跳過（.done標記）
# ============================================================
set -e
RABBIT="$1"; CHAIN="$2"
[ -z "$RABBIT" ] || [ -z "$CHAIN" ] && {
    echo "用法：bash 04d_stitch_igblast_igkl.sh <rabbit> <IGK|IGL>"; exit 1; }

BASE_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
PYTHON="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python"
ROUND_RETRY="${BASE_DIR}/script/igblast_round_retry.py"
IGBLAST_BATCH=25000

OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/${CHAIN}_Stitched"
CHUNK_DIR="${OUT_DIR}/chunks"
SAMPLE="15651-LT-${RABBIT}_AVITI_${CHAIN}_Stitched"
FINAL_AIRR="${OUT_DIR}/${SAMPLE}.airr.tsv"
STATS="${OUT_DIR}/${SAMPLE}_qc_stats.json"

GV="${IMGT}/${CHAIN}V.fasta"
GJ="${IMGT}/${CHAIN}J.fasta"
GD="${IMGT}/${CHAIN}J.fasta"

GREEN='\033[0;32m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }

split_fasta() {
    awk -v batch="$3" -v prefix="$2" '
    /^>/ { if (n % batch == 0) { if (fi > 0) close(outfile); outfile = sprintf("%s_%04d.fasta", prefix, fi); fi++ } n++ }
    { print > outfile }' "$1"
}

# 2026-09-30新增：最終合併檔已存在就直接跳過，不重新執行。
# 背景：chunk層級的中間檔案（chunks/chunk_*.airr.tsv、*.fasta）在確認最終
# 合併檔正確無誤後，可能已被手動清理以回收硬碟空間。若此時該rabbit/chain
# 又被重新呼叫這支腳本，底下「合併所有chunk」那一步沒有保護機制，會抓不到
# chunk檔案而報錯、或在最壞情況下用空/錯誤內容覆蓋掉這份已完成的正確結果。
# 加這個檢查是為了讓腳本本身具備冪等性（idempotent），可以安全地重複呼叫。
if [ -s "$FINAL_AIRR" ]; then
    info "最終合併檔已存在（${FINAL_AIRR}），視為此rabbit/chain已完成，略過重新執行。"
    info "若確定要強制重跑，請先手動搬移或刪除該檔案，並確認chunks/底下的中間檔案仍完整存在。"
    exit 0
fi

[ -f "$ROUND_RETRY" ] || error "找不到 ${ROUND_RETRY}"
[ -d "$CHUNK_DIR" ] || error "找不到 ${CHUNK_DIR}，請先執行 04c_stitch_qc_igkl.py ${RABBIT} ${CHAIN}"
shopt -s nullglob
FASTAS=("${CHUNK_DIR}"/chunk_*.fasta)
[ ${#FASTAS[@]} -gt 0 ] || error "找不到chunk fasta"
info "Rabbit ${RABBIT} ${CHAIN} Stitched：找到 ${#FASTAS[@]} 個 chunk fasta"

for FASTA in "${FASTAS[@]}"; do
    IDX=$(basename "$FASTA" .fasta | sed 's/chunk_//')
    AIRR="${CHUNK_DIR}/chunk_${IDX}.airr.tsv"
    DONE_MARK="${CHUNK_DIR}/chunk_${IDX}.igblast.done"
    SUB_DIR="${CHUNK_DIR}/chunk_${IDX}_sub"

    N_SEQ=$(grep -c "^>" "$FASTA" 2>/dev/null || echo 0)
    if [ "$N_SEQ" -eq 0 ]; then
        info "chunk ${IDX} 是空檔（沒有序列通過拼接+QC），標記完成並跳過"
        touch "$DONE_MARK"
        continue
    fi

    if [ -f "$DONE_MARK" ]; then
        info "chunk ${IDX} igblast 已完成，跳過"
        continue
    fi

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
: > "$FINAL_AIRR"
FIRST=1
for AIRR in "${AIRRS[@]}"; do
    [ -f "$AIRR" ] || continue
    if [ "$FIRST" = "1" ]; then cat "$AIRR" >> "$FINAL_AIRR"; FIRST=0; else tail -n +2 "$AIRR" >> "$FINAL_AIRR"; fi
done
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
               'method': 'coordinate-stitching + Biopython QC + round-based retry igblast',
               'reads_igblast': n, 'productive_count': prod}, f, indent=2)
print(f'productive={prod/n*100:.1f}%' if n else 'n=0')
"
info "✓ Rabbit ${RABBIT} ${CHAIN} Stitched igblast 全部完成：${FINAL_AIRR}"

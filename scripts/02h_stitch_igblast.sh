#!/bin/bash
# ============================================================
# 02h_stitch_igblast.sh
# 接續 02g_stitch_qc_igh.py：對每個拼接+QC通過的chunk fasta跑igblastn
# （完整V(D)J序列，判定跟merge/R2-only一樣完整的v_call/d_call/j_call/
# junction_aa/productive），逐塊執行、逐塊完成即算入最終檔案，最後合併
# 成一個AIRR TSV。
#
# 這是IGH的第三輪igblast（R1、R2各自先跑過一次，這裡是拼接後的序列
# 再跑一次），目的是用完整涵蓋FR1到J的序列做V基因判定，修正R2-only
# 看不到FR1/CDR1導致IGHV1S40/1S45這類近乎同源基因對判錯的問題
# （2026-09-24~25驗證，見討論記錄）。
#
# 用法：bash 02h_stitch_igblast.sh <rabbit_number>
# 長時間執行請用：
#   caffeinate -i bash 02h_stitch_igblast.sh 1 2>&1 | tee ~/Desktop/aviti_02h_stitch.log
#
# 中斷後重跑：已完成的 chunk 會被跳過（chunks/chunk_XXX.igblast.done 標記）。
# ============================================================

RABBIT="$1"
[ -z "$RABBIT" ] && { echo "用法：bash 02h_stitch_igblast.sh <rabbit_number>"; exit 1; }

BASE_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
IGDATA="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/share/igblast"
IGBLASTN="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/igblastn"
AUX="${IGDATA}/optional_file/rabbit_gl.aux"
BLASTDB="${IGDATA}/database"
THREADS=2
IGBLAST_BATCH=25000
PYTHON="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python"

OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/IGH_Stitched"
CHUNK_DIR="${OUT_DIR}/chunks"
SAMPLE="15651-LT-${RABBIT}_AVITI_IGH_Stitched"
FINAL_AIRR="${OUT_DIR}/${SAMPLE}.airr.tsv"
STATS="${OUT_DIR}/${SAMPLE}_qc_stats.json"

set -e
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }

IGBLAST_ARGS=(
    -organism rabbit -ig_seqtype Ig
    -outfmt 19 -num_threads "$THREADS"
    -extend_align5end -extend_align3end
    -auxiliary_data "$AUX"
    -germline_db_V "${IMGT}/IGHV.fasta"
    -germline_db_D "${IMGT}/IGHD.fasta"
    -germline_db_J "${IMGT}/IGHJ.fasta"
    -min_D_match 5
)

# 跟 02d/02f 完全相同的驗證＋切半重跑機制
igblast_verified() {
    local FASTA="$1" OUT="$2"
    local N_IN
    N_IN=$(grep -c "^>" "$FASTA")
    local TMP_OUT="${OUT}.attempt"

    IGDATA="$IGDATA" BLASTDB="$BLASTDB" "$IGBLASTN" "${IGBLAST_ARGS[@]}" \
        -query "$FASTA" -out "$TMP_OUT" 2>/dev/null

    local N_OUT
    N_OUT=$(tail -n +2 "$TMP_OUT" 2>/dev/null | wc -l | tr -d ' ')

    if [ "${N_OUT:-0}" -eq "$N_IN" ]; then
        mv "$TMP_OUT" "$OUT"
        return
    fi

    if [ "$N_IN" -le 1 ]; then
        warn "  無法處理的序列（igblastn持續失敗），記錄為缺失並跳過：$(grep '^>' "$FASTA")"
        mv "$TMP_OUT" "$OUT"
        return
    fi

    warn "  批次驗證失敗（輸入${N_IN}條，輸出${N_OUT:-0}條），切半重跑：$(basename "$FASTA")"
    local HALF=$(( (N_IN + 1) / 2 ))
    local SPLIT_DIR
    SPLIT_DIR=$(mktemp -d "$(dirname "$FASTA")/retry.XXXXXX")
    split_fasta "$FASTA" "${SPLIT_DIR}/part" "$HALF"

    local MERGED="${OUT}.merged"
    : > "$MERGED"
    local FIRST=1
    shopt -s nullglob
    for HALF_FASTA in "${SPLIT_DIR}"/part_*.fasta; do
        local HALF_OUT="${HALF_FASTA}.airr.tsv"
        igblast_verified "$HALF_FASTA" "$HALF_OUT"
        if [ "$FIRST" = "1" ]; then
            cat "$HALF_OUT" >> "$MERGED"
            FIRST=0
        else
            tail -n +2 "$HALF_OUT" >> "$MERGED"
        fi
    done
    rm -rf "$SPLIT_DIR"
    mv "$MERGED" "$OUT"
    rm -f "$TMP_OUT"
}

split_fasta() {
    awk -v batch="$3" -v prefix="$2" '
    /^>/ {
        if (n % batch == 0) {
            if (fi > 0) close(outfile)
            outfile = sprintf("%s_%04d.fasta", prefix, fi)
            fi++
        }
        n++
    }
    { print > outfile }
    ' "$1"
}

# 2026-09-30新增：最終合併檔已存在就直接跳過，不重新執行（避免chunk中間
# 檔案被清理後，重跑時「合併所有chunk」步驟報錯或覆蓋掉已完成的結果）。
if [ -s "$FINAL_AIRR" ]; then
    info "最終合併檔已存在（${FINAL_AIRR}），視為此rabbit已完成，略過重新執行。"
    info "若確定要強制重跑，請先手動搬移或刪除該檔案，並確認chunks/底下的中間檔案仍完整存在。"
    exit 0
fi

step "環境確認"
[ -f "$IGBLASTN" ] || error "找不到 igblastn：$IGBLASTN"
[ -d "${IGDATA}/internal_data/rabbit" ] || error "rabbit internal_data 不存在"
[ -d "$CHUNK_DIR" ] || error "找不到 ${CHUNK_DIR}，請先執行 02g_stitch_qc_igh.py ${RABBIT}"
shopt -s nullglob
FASTAS=("${CHUNK_DIR}"/chunk_*.fasta)
[ ${#FASTAS[@]} -gt 0 ] || error "找不到任何 chunk fasta，請先執行 02g_stitch_qc_igh.py ${RABBIT}"
info "找到 ${#FASTAS[@]} 個 chunk fasta"

for FASTA in "${FASTAS[@]}"; do
    IDX=$(basename "$FASTA" .fasta | sed 's/chunk_//')
    AIRR="${CHUNK_DIR}/chunk_${IDX}.airr.tsv"
    DONE_MARK="${CHUNK_DIR}/chunk_${IDX}.igblast.done"
    SUB_DIR="${CHUNK_DIR}/chunk_${IDX}_sub"

    N_SEQ=$(grep -c "^>" "$FASTA" 2>/dev/null || echo 0)
    if [ "$N_SEQ" -eq 0 ]; then
        info "chunk ${IDX} 沒有任何序列（QC全部排除），跳過"
        touch "$DONE_MARK"
        continue
    fi

    if [ -f "$DONE_MARK" ]; then
        info "chunk ${IDX} igblast 已完成，跳過"
        continue
    fi

    step "IgBLAST chunk ${IDX}（${N_SEQ} 條序列，切成 ${IGBLAST_BATCH}/批）"

    mkdir -p "$SUB_DIR"
    shopt -s nullglob
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

        if [ -f "$SUB_DONE" ]; then
            info "    ${SUB_ID} 已完成，跳過"
            continue
        fi

        SUB_N=$(grep -c "^>" "$SUB_FASTA")
        info "    跑 ${SUB_ID}（${SUB_N} 條）..."

        igblast_verified "$SUB_FASTA" "$SUB_AIRR"
        SUB_OUT_N=$(tail -n +2 "$SUB_AIRR" | wc -l | tr -d ' ')
        [ "$SUB_OUT_N" -eq "$SUB_N" ] || warn "    ${SUB_ID} 驗證後仍有缺失：輸入${SUB_N}條，輸出${SUB_OUT_N}條"

        touch "$SUB_DONE"
    done

    SUB_AIRRS=("${SUB_DIR}"/sub_*.airr.tsv)
    head -1 "${SUB_AIRRS[0]}" > "$AIRR"
    for SUB_AIRR in "${SUB_AIRRS[@]}"; do
        tail -n +2 "$SUB_AIRR" >> "$AIRR"
    done
    rm -rf "$SUB_DIR"

    touch "$DONE_MARK"
    info "  ✓ chunk ${IDX} 完成：$(tail -n +2 "$AIRR" | wc -l | tr -d ' ') rows"
done

step "合併所有 chunk 的 AIRR TSV"
AIRRS=("${CHUNK_DIR}"/chunk_*.airr.tsv)
if [ ${#AIRRS[@]} -eq 0 ]; then
    error "沒有任何chunk產生AIRR TSV結果（可能所有chunk都被QC排除光了），請檢查02g的輸出"
fi
head -1 "${AIRRS[0]}" > "$FINAL_AIRR"
for AIRR in "${AIRRS[@]}"; do
    tail -n +2 "$AIRR" >> "$FINAL_AIRR"
done
TOTAL_N=$(tail -n +2 "$FINAL_AIRR" | wc -l | tr -d ' ')
info "合併完成：${FINAL_AIRR}（${TOTAL_N} rows）"

"$PYTHON" -c "
import json, csv
n = 0; prod = 0; has_j = 0
with open('${FINAL_AIRR}') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        n += 1
        if row['productive'] == 'T': prod += 1
        if row['junction_aa']: has_j += 1
with open('${STATS}', 'w') as f:
    json.dump({'sample': '${SAMPLE}', 'platform': 'AVITI', 'chain': 'IGH',
               'method': 'Coordinate-guided Stitching (R1+R2 stitched by germline coordinates) + QC (Biopython local alignment, identity>=90%)',
               'reads_igblast': n, 'productive_count': prod, 'has_junction_aa': has_j}, f, indent=2)
print(f'productive={prod/n*100:.1f}%  has_junction_aa={has_j/n*100:.1f}%')
"
info "✓ Rabbit ${RABBIT} IGH 座標拼接+QC 全部完成：${FINAL_AIRR}"

#!/bin/bash
# ============================================================
# 02e_run_aviti_r1only_igh.sh
# AVITI pipeline — IGH R1 側路徑：跟 02c（R2-only）配套，抽出同一批
# demux 命中的 read pairs 裡「R1角色」那端（TSO-anchored side，涵蓋
# FR1起～多半到CDR1/FR2一帶），供之後用read ID跟R2-only結果配對，
# 組合出比純R2-only更完整的productive判定（R2負責D/J/junction/frame，
# R1額外檢查FR1/CDR1有沒有stop codon）。
#
# 背景與完整討論見 discussion_notes/aviti_r2only_preanalysis_20260915/
# 與 2026-09-22~23 對話記錄。R1單獨幾乎讀不到CDR3（junction判讀率僅
# 0.2%），所以R1只拿來做productive的補充驗證，不取代R2作為CDR3/D/J
# 主要來源。
#
# 跟 02c_run_aviti_r2only_igh.sh 的差異：
#   - demux/rescue邏輯完全相同（同一批raw reads、同一個chunk切法，
#     確保跟02c產出的chunk一一對應、read ID可以互相配對）
#   - 保留的是「R1角色」的read（s1_R1 + s2_R2），不是R2角色
#   - R1只做TSO trim（-g T7_TSO），不做R2那種quality-cutoff/長度篩選
#     ——避免無謂地把R1的yield拉低，且經10K read驗證過此作法可行
#     （TSO trim後100%的read都保留，無因trim被拋棄的情況）
#
# 用法：bash 02e_run_aviti_r1only_igh.sh <rabbit_number>
#   caffeinate -i bash 02e_run_aviti_r1only_igh.sh 1 2>&1 | tee ~/Desktop/aviti_02e_r1.log
#
# 中斷後重跑：已完成的 chunk 會被跳過（chunks/chunk_XXX.done 標記）。
# ============================================================

RABBIT="$1"
[ -z "$RABBIT" ] && { echo "用法：bash 02e_run_aviti_r1only_igh.sh <rabbit_number>"; exit 1; }

BASE_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
R1="${BASE_DIR}/samples/15651-LT-${RABBIT}_R1.fastq"
R2="${BASE_DIR}/samples/15651-LT-${RABBIT}_R2.fastq"
OUT_DIR="${BASE_DIR}/results/15651-LT-${RABBIT}_Aviti24/IGH_R1only"
CHUNK_DIR="${OUT_DIR}/chunks"
TMP_DIR="${OUT_DIR}/tmp"
PYTHON="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python"

CHUNK_RAW_PAIRS=1500000   # 跟02c完全相同，確保chunk切法一致、read ID可對應

T7_TSO="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"
IGH_REV1="CAGTGGGAAGACTGACGGAGCCTTAG"
IGH_REV2="CAGTGGGAAGACTGATGGAGCCTTAG"

set -e
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
count_reads() { echo $(( $(wc -l < "$1") / 4 )); }

step "環境確認"
[ -f "$R1" ] || error "找不到 R1：$R1"
[ -f "$R2" ] || error "找不到 R2：$R2"
command -v cutadapt &>/dev/null || error "找不到 cutadapt（記得先 conda activate igblast-env）"
mkdir -p "$CHUNK_DIR" "$TMP_DIR"
[ -f "${OUT_DIR}/chunk_stats.csv" ] || echo "chunk,raw_pairs,igh_demux,tso_trimmed" > "${OUT_DIR}/chunk_stats.csv"

TOTAL_LINES=$(wc -l < "$R1")
TOTAL_PAIRS=$(( TOTAL_LINES / 4 ))
NUM_CHUNKS=$(( (TOTAL_PAIRS + CHUNK_RAW_PAIRS - 1) / CHUNK_RAW_PAIRS ))
info "Rabbit ${RABBIT}：raw total ${TOTAL_PAIRS} pairs，切成 ${NUM_CHUNKS} 塊（每塊上限 ${CHUNK_RAW_PAIRS} pairs）"

for (( i=0; i<NUM_CHUNKS; i++ )); do
    IDX=$(printf "%03d" "$i")
    DONE_MARK="${CHUNK_DIR}/chunk_${IDX}.done"
    OUT_FASTA="${CHUNK_DIR}/chunk_${IDX}.fasta"

    if [ -f "$DONE_MARK" ]; then
        info "chunk ${IDX} 已完成，跳過"
        continue
    fi

    step "Chunk ${IDX} / $(printf '%03d' $((NUM_CHUNKS-1)))"
    START_LINE=$(( i * CHUNK_RAW_PAIRS * 4 + 1 ))
    END_LINE=$(( (i+1) * CHUNK_RAW_PAIRS * 4 ))

    sed -n "${START_LINE},${END_LINE}p" "$R1" > "${TMP_DIR}/raw_R1.fastq"
    sed -n "${START_LINE},${END_LINE}p" "$R2" > "${TMP_DIR}/raw_R2.fastq"
    N_RAW=$(count_reads "${TMP_DIR}/raw_R1.fastq")
    if [ "$N_RAW" -eq 0 ]; then
        info "chunk ${IDX} 沒有資料（已到檔案結尾），結束"
        rm -f "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"
        break
    fi
    info "  raw: ${N_RAW} pairs"

    # Step1：demux（primer直接命中在raw R2上，跟02c完全相同邏輯）
    cutadapt -G "$IGH_REV1" -G "$IGH_REV2" -e 0.1 --overlap 15 \
        -o "${TMP_DIR}/s1_R1.fastq" -p "${TMP_DIR}/s1_R2.fastq" \
        --untrimmed-output "${TMP_DIR}/non_R1.fastq" \
        --untrimmed-paired-output "${TMP_DIR}/non_R2.fastq" \
        --quiet "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"
    rm -f "${TMP_DIR}/raw_R1.fastq" "${TMP_DIR}/raw_R2.fastq"

    # Step2：rescue（primer落在raw R1上，方向對調）
    cutadapt -g "$IGH_REV1" -g "$IGH_REV2" -e 0.1 --overlap 15 \
        -o "${TMP_DIR}/s2_trimR1.fastq" -p "${TMP_DIR}/s2_R2.fastq" \
        --untrimmed-output /dev/null --untrimmed-paired-output /dev/null \
        --quiet "${TMP_DIR}/non_R1.fastq" "${TMP_DIR}/non_R2.fastq"
    rm -f "${TMP_DIR}/non_R1.fastq" "${TMP_DIR}/non_R2.fastq"

    # Step3：彙整兩個來源的「R1角色」read。
    # s1組：primer命中在s1_R2，所以s1_R1本來就是正確方向的R1角色。
    # s2組（rescue，方向對調）：primer命中在s2_trimR1（=R2角色），
    # 所以s2_R2才是R1角色（對應02c那邊拿s2_trimR1當R2角色是同一組對調邏輯）。
    # 這樣組出來的combined_R1，跟02c的combined_R2是完全對應的同一批read pairs。
    cat "${TMP_DIR}/s1_R1.fastq" "${TMP_DIR}/s2_R2.fastq" > "${TMP_DIR}/combined_R1.fastq"
    rm -f "${TMP_DIR}/s1_R1.fastq" "${TMP_DIR}/s1_R2.fastq" "${TMP_DIR}/s2_R2.fastq" "${TMP_DIR}/s2_trimR1.fastq"
    N_IGH=$(count_reads "${TMP_DIR}/combined_R1.fastq")
    info "  IGH demux: ${N_IGH} reads"

    # R1 TSO trim（不做quality-cutoff/長度篩選，經10K read驗證此作法
    # 不會無謂丟棄read；序列品質交給igblast的v_identity/比對結果判斷）
    cutadapt -g "$T7_TSO" -e 0.1 --overlap 15 --quiet \
        "${TMP_DIR}/combined_R1.fastq" -o "${TMP_DIR}/trimmed.fastq"
    rm -f "${TMP_DIR}/combined_R1.fastq"
    N_TRIM=$(count_reads "${TMP_DIR}/trimmed.fastq")
    info "  TSO trim 後: ${N_TRIM} reads"

    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "${TMP_DIR}/trimmed.fastq" > "$OUT_FASTA"
    rm -f "${TMP_DIR}/trimmed.fastq"

    echo "${IDX},${N_RAW},${N_IGH},${N_TRIM}" >> "${OUT_DIR}/chunk_stats.csv"
    touch "$DONE_MARK"
    info "  ✓ chunk ${IDX} 完成：${N_RAW} raw → ${N_IGH} IGH demux → ${N_TRIM} TSO trim 後"
done

rmdir "$TMP_DIR" 2>/dev/null || true
step "Demux + R1-only TSO trim 完成"
TOTAL_FASTA=$(cat "${CHUNK_DIR}"/chunk_*.fasta 2>/dev/null | grep -c "^>")
N_CHUNKS_DONE=$(ls "${CHUNK_DIR}"/chunk_*.fasta 2>/dev/null | wc -l | tr -d ' ')
info "Rabbit ${RABBIT} IGH R1-only 總計：${TOTAL_FASTA} 條序列，${N_CHUNKS_DONE} 個 chunk fasta"
info "接下來執行：bash 02f_run_aviti_r1only_igblast.sh ${RABBIT}"

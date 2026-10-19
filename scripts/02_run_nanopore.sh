#!/bin/bash
IGDATA="${BCR_IGDATA:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/share/igblast}"
PYTHON="${BCR_PIPELINE_PYTHON:-/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python}"
BASE_DIR="${BCR_PIPELINE_ROOT:-/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline}"
SAMPLES_DIR="${BASE_DIR}/samples/MZPD8T_Nanopore"
IMGT="${BASE_DIR}/reference/IMGT_reference_clean"
WORK_DIR="${BASE_DIR}/results/MZPD8T_Nanopore"
OUT_DIR="${WORK_DIR}/results"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLASTDB="${IGDATA}/database"
AUX="${IGDATA}/optional_file/rabbit_gl.aux"
IGBLASTN="${IGDATA}/bin/igblastn"
THREADS=2

# WT 樣本
WT_IGH_FASTQ="${SAMPLES_DIR}/MZPD8T_1_WT-IGHV.fastq"
WT_IGK_FASTQ="${SAMPLES_DIR}/MZPD8T_2_WT-IGK.fastq"
WT_IGL_FASTQ="${SAMPLES_DIR}/MZPD8T_3_WT-IGL.fastq"

# Immu 樣本
IMMU_IGH_FASTQ="${SAMPLES_DIR}/MZPD8T_4_Immu-IGHV.fastq"
IMMU_IGK_FASTQ="${SAMPLES_DIR}/MZPD8T_5_Immu-IGK.fastq"
IMMU_IGL_FASTQ="${SAMPLES_DIR}/MZPD8T_6_Immu-IGL.fastq"

TSO="AAGCAGTGGTATCAACGCAGAGT"
IGH_REV1="CAGTGGGAAGACTGACGGAGCCTTAG"
IGH_REV2="CAGTGGGAAGACTGATGGAGCCTTAG"
IGK_REV1="TGGTGGGAAGAKGAGGACAGTAGG"
IGK_REV2="TGGTGGGAAGAKGAGGACACTAGG"
IGK_REV3="TGGTGGGAAGAKGAGGACAGAAGG"
IGL_REV1="CAAGGGGGCGACCACAGGCTGAC"
IGL_REV2="GTGAAGGAGTGACTACGGGTTGACC"
IGL_REV3="GAGGGGGTCACCGCGGGCTGAC"

LEN_MIN_IGH=350; LEN_MIN_IGK=320; LEN_MIN_IGL=320; LEN_MAX=600

set -e
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; CYAN='\033[0;36m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step()  { echo -e "\n${CYAN}━━━ $1 ━━━${NC}"; }
count_reads() { echo $(( $(wc -l < "$1") / 4 )); }
rc_iupac() {
    "$PYTHON" -c "
s = '$1'.upper()
comp = {'A':'T','T':'A','G':'C','C':'G','N':'N','K':'M','M':'K','R':'Y','Y':'R','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D'}
print(''.join(comp.get(b,'N') for b in reversed(s)))
"
}

step "環境確認"
command -v chopper &>/dev/null || error "找不到 chopper"
command -v cutadapt &>/dev/null || error "找不到 cutadapt"
[ -f "$IGBLASTN" ]                      || error "找不到 igblastn：$IGBLASTN"
[ -d "${IGDATA}/internal_data/rabbit" ] || error "rabbit internal_data 不存在"
[ -f "$AUX" ] && info "rabbit_gl.aux: ✓" || warn "rabbit_gl.aux 未找到"
mkdir -p "$OUT_DIR"
info "IGDATA：$IGDATA"
info "igblastn：$($IGBLASTN -version 2>&1 | head -1)"

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

# igblastn 1.22.0 已知問題：單次呼叫累積處理夠多序列後會對特定序列丟出
# CEXCEPTION(NULL pointer)，該序列之後排隊的部分被靜默丟棄且exit code仍是0
# （2026-09-17 於 AVITI IGK/IGL 實測發現）。這裡核對輸入輸出行數，對不上就切半
# 重跑，保證不再靜默漏資料。Nanopore單一樣本讀數不多，通常一次就會過。
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

run_nano_chain() {
    local CHAIN="$1" LOCUS="$2" FASTQ="$3" LEN_MIN="$4" GROUP="$5"
    shift 5
    local REV_PRIMERS=("$@")
    local SAMPLE="MZPD8T_${GROUP}_${CHAIN}"
    local SDIR="${OUT_DIR}/${SAMPLE}"
    mkdir -p "$SDIR"

    step "處理 ${SAMPLE}"
    if [ -f "${SDIR}/${SAMPLE}_igblast.airr.tsv" ]; then
        info "${SAMPLE} 先前已完成，跳過"
        return
    fi
    [ -f "$FASTQ" ] || error "找不到：$FASTQ"
    local RAW; RAW=$(count_reads "$FASTQ")
    info "  Raw reads：${RAW}"

    info "Step 1/4 — chopper (-q 10, ${LEN_MIN}–${LEN_MAX} bp)"
    local QC="${SDIR}/${CHAIN}_qc.fastq"
    chopper -q 10 -l "$LEN_MIN" --maxlength "$LEN_MAX" < "$FASTQ" > "$QC"
    local QC_N; QC_N=$(count_reads "$QC")
    info "  After QC：${QC_N} ($("$PYTHON" -c "print(f'{${QC_N}/${RAW}*100:.1f}%')"))"

    info "Step 2/4 — Adapter Trim (-e 0.1, --revcomp)"
    local TRIM="${SDIR}/${CHAIN}_trimmed.fastq"
    local CA_ARGS=(-g "$TSO")
    for REV in "${REV_PRIMERS[@]}"; do
        CA_ARGS+=(-a "$(rc_iupac "$REV")")
    done
    CA_ARGS+=(--revcomp -e 0.1 --overlap 15 --discard-untrimmed
              --minimum-length "$LEN_MIN" --maximum-length "$LEN_MAX" --quiet)
    cutadapt "${CA_ARGS[@]}" "$QC" -o "$TRIM"
    local TRIM_N; TRIM_N=$(count_reads "$TRIM")
    info "  After trim：${TRIM_N} ($("$PYTHON" -c "print(f'{${TRIM_N}/${QC_N}*100:.1f}%')"))"

    info "Step 3/4 — FASTQ → FASTA"
    local FASTA="${SDIR}/${CHAIN}_trimmed.fasta"
    awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$TRIM" > "$FASTA"

    info "Step 4/4 — IgBLAST"
    local AIRR="${SDIR}/${SAMPLE}_igblast.airr.tsv"
    local IGBLAST_ARGS=(
        -organism rabbit -ig_seqtype Ig
        -outfmt 19 -num_threads "$THREADS"
        -extend_align5end -extend_align3end
        -auxiliary_data "$AUX"
    )

    if [ "$LOCUS" = "IGH" ]; then
        IGBLAST_ARGS+=(
            -germline_db_V "${IMGT}/IGHV.fasta"
            -germline_db_D "${IMGT}/IGHD.fasta"
            -germline_db_J "${IMGT}/IGHJ.fasta"
            -min_D_match 5)
    elif [ "$LOCUS" = "IGK" ]; then
        IGBLAST_ARGS+=(
            -germline_db_V "${IMGT}/IGKV.fasta"
            -germline_db_D "${IMGT}/IGKJ.fasta"
            -germline_db_J "${IMGT}/IGKJ.fasta")
    elif [ "$LOCUS" = "IGL" ]; then
        IGBLAST_ARGS+=(
            -germline_db_V "${IMGT}/IGLV.fasta"
            -germline_db_D "${IMGT}/IGLJ.fasta"
            -germline_db_J "${IMGT}/IGLJ.fasta")
    fi

    igblast_verified "$FASTA" "$AIRR"

    local IGBLAST_N
    IGBLAST_N=$(tail -n +2 "$AIRR" | wc -l | tr -d ' ')
    [ "$IGBLAST_N" -eq "$TRIM_N" ] || warn "  驗證後仍有缺失：輸入${TRIM_N}條，輸出${IGBLAST_N}條"
    info "  IgBLAST 輸出：${IGBLAST_N} sequences"

    "$PYTHON" -c "
import json
with open('${SDIR}/${SAMPLE}_qc_stats.json','w') as f:
    json.dump({'sample':'${SAMPLE}','platform':'Nanopore','chain':'${LOCUS}',
               'reads_raw':${RAW},'reads_qc':${QC_N},
               'reads_trimmed':${TRIM_N},'reads_igblast':${IGBLAST_N}}, f, indent=2)
"
    info "✓ ${SAMPLE} 完成"
}

run_nano_chain "IGHV" "IGH" "$WT_IGH_FASTQ"   "$LEN_MIN_IGH" "WT"   "$IGH_REV1" "$IGH_REV2"
run_nano_chain "IGK"  "IGK" "$WT_IGK_FASTQ"   "$LEN_MIN_IGK" "WT"   "$IGK_REV1" "$IGK_REV2" "$IGK_REV3"
run_nano_chain "IGL"  "IGL" "$WT_IGL_FASTQ"   "$LEN_MIN_IGL" "WT"   "$IGL_REV1" "$IGL_REV2" "$IGL_REV3"
run_nano_chain "IGHV" "IGH" "$IMMU_IGH_FASTQ" "$LEN_MIN_IGH" "Immu" "$IGH_REV1" "$IGH_REV2"
run_nano_chain "IGK"  "IGK" "$IMMU_IGK_FASTQ" "$LEN_MIN_IGK" "Immu" "$IGK_REV1" "$IGK_REV2" "$IGK_REV3"
run_nano_chain "IGL"  "IGL" "$IMMU_IGL_FASTQ" "$LEN_MIN_IGL" "Immu" "$IGL_REV1" "$IGL_REV2" "$IGL_REV3"

echo -e "\n${GREEN}============================================${NC}"
echo -e "${GREEN}  Nanopore pipeline 完成！輸出：${OUT_DIR}${NC}"
echo -e "${GREEN}============================================${NC}"
ls -lh "${OUT_DIR}"/MZPD8T_*/*.airr.tsv 2>/dev/null || true

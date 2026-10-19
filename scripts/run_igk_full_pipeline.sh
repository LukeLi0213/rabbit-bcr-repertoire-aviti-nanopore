#!/bin/bash
# IGK全量座標拼接+QC pipeline正式執行入口（只跑IGK，IGL暫停）。
# 3隻兔子，每隻依序：demux(04a) -> R1-only igblast(04b) -> R2-only igblast(04b)
# -> 座標拼接+QC(04c) -> 第3次igblast(04d)。
# 中斷後重跑：各步驟自身都有chunk/子批次層級的.done標記，會自動跳過已完成部分。
set -e
SCRIPT_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline/script"
PYTHON="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python"
CHAIN="IGK"

for RABBIT in 1 2 3; do
    echo ""
    echo "################################################################"
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step1 demux 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    echo "################################################################"
    bash "${SCRIPT_DIR}/04a_run_aviti_igkl_demux.sh" "$RABBIT" "$CHAIN"

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step2 R1-only igblast 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    bash "${SCRIPT_DIR}/04b_run_aviti_igkl_igblast.sh" "$RABBIT" "$CHAIN" r1only

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step2 R2-only igblast 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    bash "${SCRIPT_DIR}/04b_run_aviti_igkl_igblast.sh" "$RABBIT" "$CHAIN" r2only

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step3 座標拼接+QC 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    "$PYTHON" "${SCRIPT_DIR}/04c_stitch_qc_igkl.py" "$RABBIT" "$CHAIN"

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — Step4 拼接後第3次igblast 開始 $(date '+%Y-%m-%d %H:%M:%S')"
    bash "${SCRIPT_DIR}/04d_stitch_igblast_igkl.sh" "$RABBIT" "$CHAIN"

    echo "### Rabbit ${RABBIT} ${CHAIN} 全部完成 $(date '+%Y-%m-%d %H:%M:%S')"
done

echo ""
echo "IGK_FULL_PIPELINE_ALL_DONE"

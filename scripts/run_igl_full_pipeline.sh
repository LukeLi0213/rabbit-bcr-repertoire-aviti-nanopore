#!/bin/bash
# IGL全量座標拼接+QC pipeline正式執行入口。
# 3隻兔子，每隻依序：demux(04a) -> R1-only igblast(04b) -> R2-only igblast(04b)
# -> 座標拼接+QC(04c) -> 第3次igblast(04d)。
#
# 已知背景（見記憶bcr_pipeline_igl_r2only_retry_slowdown.md）：IGL R2-only/
# 拼接後第3次igblast會頻繁觸發igblastn 1.22.0本身的NULL pointer bug（已確認
# 是igblastn的問題，不是我們資料或pipeline的問題，查不到不犧牲準確度的修正
# 方式），round-based capped retry會自動處理、但速度慢很多，全量預估6-9天。
#
# 2026-10-02新增：跟IGK當初做法不同，這次**每隻兔子跑完Step4後，確認3份
# 最終合併檔存在且非空，就立即刪除該兔子的chunks/資料夾**，不等三隻兔子
# 都跑完才清——因為IGL全程耗時遠比IGK長，硬碟空間需要全程控制在「只疊加
# 當下這一隻」的量，不要累積到三隻份。
#
# 中斷後重跑：各步驟自身都有chunk/子批次層級的.done標記，會自動跳過已完成
# 部分；但注意如果某隻兔子的chunks已經被清掉、該兔子理論上不會再被重跑
# （04b/04d已有「最終檔案存在就跳過」的保護），所以中斷後重跑是安全的。
set -e
SCRIPT_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline/script"
RESULTS_DIR="/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline/results"
PYTHON="/Users/newuser/Desktop/claude-code/0-Environment/miniforge3/envs/igblast-env/bin/python"
CHAIN="IGL"

cleanup_rabbit_chunks() {
    local rabbit="$1"
    local base="${RESULTS_DIR}/15651-LT-${rabbit}_Aviti24"
    local ok=1

    for side_dir in "${CHAIN}_R1only" "${CHAIN}_R2only" "${CHAIN}_Stitched"; do
        local final_airr
        final_airr=$(ls "${base}/${side_dir}"/*.airr.tsv 2>/dev/null | head -1)
        if [ -z "$final_airr" ] || [ ! -s "$final_airr" ]; then
            echo "  !! 警告：${side_dir} 的最終合併檔不存在或是空檔，不會清理 chunks/（${base}/${side_dir}）"
            ok=0
        fi
    done

    if [ "$ok" -eq 1 ]; then
        for side_dir in "${CHAIN}_R1only" "${CHAIN}_R2only" "${CHAIN}_Stitched"; do
            local chunks_dir="${base}/${side_dir}/chunks"
            if [ -d "$chunks_dir" ]; then
                local size
                size=$(du -sh "$chunks_dir" 2>/dev/null | cut -f1)
                rm -rf "$chunks_dir"
                echo "  已清理: ${side_dir}/chunks (${size})"
            fi
        done
    else
        echo "  !! Rabbit ${rabbit} 的chunks/保留不清，請人工確認後再手動處理"
    fi
}

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

    echo ""
    echo "### Rabbit ${RABBIT} ${CHAIN} — 驗證最終檔案並清理chunks $(date '+%Y-%m-%d %H:%M:%S')"
    cleanup_rabbit_chunks "$RABBIT"

    echo "### Rabbit ${RABBIT} ${CHAIN} 全部完成 $(date '+%Y-%m-%d %H:%M:%S')"
done

echo ""
echo "IGL_FULL_PIPELINE_ALL_DONE"

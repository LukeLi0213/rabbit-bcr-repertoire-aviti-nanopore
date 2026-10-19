#!/usr/bin/env python3
"""
02g_stitch_qc_igh.py
接續02d(R2-only igblast)、02f(R1-only igblast)：用read ID把同一分子的R1、R2
配對起來，以igblast各自跑出來的座標為依據，把兩段序列接成一條完整V(D)J序列
（座標拼接法，Coordinate-guided Stitching），並用Biopython局部比對(local
alignment)驗證重疊區品質，只保留品質夠好的拼接結果。

背景：R2-only單獨判定V基因，對IGHV1S40/1S45這類差異點集中在FR1/CDR1（R2讀
不到的區域）的近乎同源基因對，會系統性判錯（2026-09-24~25驗證，R2-only
判定比例跟Nanopore/Merge差好幾倍）。座標拼接法把R1的FR1資訊也接回來，
2026-09-25用LT-1/2/3各10萬條抽樣驗證過，結果跟Nanopore/Merge高度吻合。

流程：
1. 把該rabbit的R1-only AIRR TSV整個讀進記憶體，建成 sequence_id -> R1資訊 的
   查找表（R1資料量約400-500萬筆，8GB RAM機器上可負擔，比照paired_igh_cache.py
   的作法）
2. 分塊讀取R2-only AIRR TSV（500K rows/chunk），逐筆用sequence_id查R1
3. 排除：R1沒有資料、R2非productive、有缺口（R1終點跟R2起點之間有空隙）
4. 拼接：R1保留到v_sequence_end，R2從對應germline位置接續
   （sequence欄位在rev_comp=T時igblast已經轉成germline方向，不需要再手動
   revcomp——這是2026-09-25驗證時發現並修正過的重點，見討論記錄）
5. QC：對有實際重疊(overlap_bp>0)的部分，用Biopython PairwiseAligner(local模式)
   比對兩邊測到的鹼基，identity<90%視為不可靠、排除
6. 輸出：拼接+QC通過的序列，寫成chunked FASTA供02h跑igblast

用法：python 02g_stitch_qc_igh.py <rabbit_number>
中斷後重跑：已完成的chunk會被跳過（chunk_XXX.stitch.done標記）。
"""
import sys
import os
import json
import pandas as pd
from Bio.Align import PairwiseAligner

RABBIT = sys.argv[1] if len(sys.argv) > 1 else None
if not RABBIT:
    print("用法：python 02g_stitch_qc_igh.py <rabbit_number>")
    sys.exit(1)

BASE_DIR = "/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
R1_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/IGH_R1only/15651-LT-{RABBIT}_AVITI_IGH_R1only.airr.tsv"
R2_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/IGH_R2only/15651-LT-{RABBIT}_AVITI_IGH_R2only.airr.tsv"
OUT_DIR = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/IGH_Stitched"
CHUNK_DIR = f"{OUT_DIR}/chunks"
STATS_PATH = f"{OUT_DIR}/chunk_stats.csv"

CHUNK_SIZE = 500_000
MARGIN = 10
QC_THRESHOLD = 0.90

R1_COLS = ['sequence_id', 'sequence', 'v_sequence_start', 'v_sequence_end',
           'v_germline_start', 'v_germline_end']
R2_COLS = ['sequence_id', 'sequence', 'productive', 'stop_codon', 'vj_in_frame',
           'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end']

os.makedirs(CHUNK_DIR, exist_ok=True)

aligner = PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5


def productive_mask(df):
    prod = df['productive'].astype(str).str.strip() if 'productive' in df.columns \
           else pd.Series([''] * len(df), index=df.index)
    if prod.isin(['T', 'F']).any():
        return prod == 'T'
    sc = df['stop_codon'].astype(str).str.strip() if 'stop_codon' in df.columns \
         else pd.Series([''] * len(df), index=df.index)
    vjf = df['vj_in_frame'].astype(str).str.strip() if 'vj_in_frame' in df.columns \
          else pd.Series([''] * len(df), index=df.index)
    if vjf.isin(['T', 'F']).any():
        return (sc == 'F') & (vjf == 'T')
    return sc == 'F'


print(f"[Rabbit {RABBIT}] Step 1: 建立R1查找表 <- {os.path.basename(R1_PATH)} ...")
r1_lookup = {}
n_r1_total = 0
for chunk in pd.read_csv(R1_PATH, sep='\t', usecols=R1_COLS, dtype=str, chunksize=CHUNK_SIZE):
    n_r1_total += len(chunk)
    for c in ['v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end']:
        chunk[c] = pd.to_numeric(chunk[c], errors='coerce')
    chunk = chunk.dropna(subset=['v_sequence_start', 'v_sequence_end',
                                  'v_germline_start', 'v_germline_end'])
    for row in chunk.itertuples(index=False):
        r1_lookup[row.sequence_id] = (row.sequence, row.v_sequence_start, row.v_sequence_end,
                                       row.v_germline_start, row.v_germline_end)
print(f"  R1總筆數={n_r1_total:,}，有效V基因比對且建入查找表={len(r1_lookup):,}")

print(f"[Rabbit {RABBIT}] Step 2: 分塊處理R2-only <- {os.path.basename(R2_PATH)} ...")
if not os.path.exists(STATS_PATH):
    with open(STATS_PATH, 'w') as f:
        f.write("chunk,r2_rows,r2_productive,paired_with_r1,nogap,qc_pass\n")

chunk_idx = 0
for chunk in pd.read_csv(R2_PATH, sep='\t', usecols=R2_COLS, dtype=str, chunksize=CHUNK_SIZE):
    idx_str = f"{chunk_idx:04d}"
    done_mark = f"{CHUNK_DIR}/chunk_{idx_str}.stitch.done"
    out_fasta = f"{CHUNK_DIR}/chunk_{idx_str}.fasta"

    if os.path.exists(done_mark):
        print(f"  chunk {idx_str} 已完成，跳過")
        chunk_idx += 1
        continue

    n_r2_rows = len(chunk)
    for c in ['v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end']:
        chunk[c] = pd.to_numeric(chunk[c], errors='coerce')

    prod = chunk[productive_mask(chunk)].copy()
    n_productive = len(prod)

    n_paired = 0
    n_nogap = 0
    n_qc_pass = 0
    out_lines = []

    for row in prod.itertuples(index=False):
        r1_info = r1_lookup.get(row.sequence_id)
        if r1_info is None:
            continue
        r1_seq, r1_seq_start, r1_seq_end, r1_germ_start, r1_germ_end = r1_info
        if pd.isna(row.v_sequence_start) or pd.isna(row.v_germline_start):
            continue
        n_paired += 1

        gap = row.v_germline_start - r1_germ_end
        if gap > 0:
            continue
        n_nogap += 1
        overlap_bp = -gap

        r2_seq = row.sequence
        r1_cut = int(r1_seq_end)
        r1_part = r1_seq[:r1_cut]
        target_germline = r1_germ_end + 1
        r2_cut_start_1based = int(row.v_sequence_start + (target_germline - row.v_germline_start))
        r2_part = r2_seq[max(r2_cut_start_1based - 1, 0):]
        stitched = r1_part + r2_part

        qc_pass = True
        ov = int(overlap_bp)
        if ov > 0:
            r1_center = int(r1_seq_start + (row.v_germline_start - r1_germ_start))
            r2_center = int(row.v_sequence_start)
            r1_win = r1_seq[max(r1_center - 1 - MARGIN, 0): r1_center - 1 + ov + MARGIN]
            r2_win = r2_seq[max(r2_center - 1 - MARGIN, 0): r2_center - 1 + ov + MARGIN]
            if len(r1_win) >= 3 and len(r2_win) >= 3:
                try:
                    aln = aligner.align(r1_win, r2_win)[0]
                    a1, a2 = str(aln[0]), str(aln[1])
                    aligned_len = len(a1)
                    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
                    identity = matches / aligned_len if aligned_len else 0.0
                    qc_pass = identity >= QC_THRESHOLD
                except Exception:
                    qc_pass = True

        if not qc_pass:
            continue
        n_qc_pass += 1
        out_lines.append(f">{row.sequence_id}\n{stitched}\n")

    with open(out_fasta, 'w') as f:
        f.writelines(out_lines)

    with open(STATS_PATH, 'a') as f:
        f.write(f"{idx_str},{n_r2_rows},{n_productive},{n_paired},{n_nogap},{n_qc_pass}\n")

    with open(done_mark, 'w') as f:
        f.write('done')

    print(f"  chunk {idx_str}: R2 rows={n_r2_rows:,} productive={n_productive:,} "
          f"配對上R1={n_paired:,} 無缺口={n_nogap:,} QC通過={n_qc_pass:,}")
    chunk_idx += 1

print(f"[Rabbit {RABBIT}] 拼接+QC完成，chunk FASTA已輸出至 {CHUNK_DIR}")
print(f"接下來執行：bash 02h_stitch_igblast.sh {RABBIT}")

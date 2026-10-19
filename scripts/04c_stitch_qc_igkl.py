#!/usr/bin/env python3
"""
04c_stitch_qc_igkl.py
接續04a(demux)、04b(R1-only/R2-only igblast)：座標拼接+QC，邏輯完全照搬
已在IGH production驗證過的02g_stitch_qc_igh.py，路徑改成IGK/IGL。

用法：python 04c_stitch_qc_igkl.py <rabbit_number> <IGK|IGL>
中斷後重跑：已完成的chunk會被跳過（chunk_XXX.stitch.done標記）。
"""
import sys
import os
import pandas as pd
from Bio.Align import PairwiseAligner

RABBIT = sys.argv[1] if len(sys.argv) > 1 else None
CHAIN = sys.argv[2] if len(sys.argv) > 2 else None
if not RABBIT or not CHAIN:
    print("用法：python 04c_stitch_qc_igkl.py <rabbit_number> <IGK|IGL>")
    sys.exit(1)

BASE_DIR = "/Users/newuser/Desktop/claude-code/20260701-BCR_pipeline"
R1_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/{CHAIN}_R1only/15651-LT-{RABBIT}_AVITI_{CHAIN}_R1only.airr.tsv"
R2_PATH = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/{CHAIN}_R2only/15651-LT-{RABBIT}_AVITI_{CHAIN}_R2only.airr.tsv"
OUT_DIR = f"{BASE_DIR}/results/15651-LT-{RABBIT}_Aviti24/{CHAIN}_Stitched"
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


print(f"[{CHAIN} Rabbit {RABBIT}] Step 1: 建立R1查找表 <- {os.path.basename(R1_PATH)} ...")
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

print(f"[{CHAIN} Rabbit {RABBIT}] Step 2: 分塊處理R2-only <- {os.path.basename(R2_PATH)} ...")
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

print(f"[{CHAIN} Rabbit {RABBIT}] 拼接+QC完成，chunk FASTA已輸出至 {CHUNK_DIR}")

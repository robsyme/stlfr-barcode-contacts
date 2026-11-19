#!/usr/bin/env python3
"""Barcode co-occurrence matrix for one chromosome.

stdin: whitespace-separated integers `A B C bin`, one line per alignment record,
where A_B_C is the stLFR barcode and bin = int(POS / bin_size) (1-based SAM POS,
as in BGI stLFRsv Stat_share.pl).

Output: `binA_start<TAB>binB_start<TAB>count`, symmetric (both orientations,
diagonal once). count = number of distinct barcodes seen in both bins.
"""
import argparse
import gzip
import json
import sys

import numpy as np
import pandas as pd
from scipy import sparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin-size", type=int, required=True)
    ap.add_argument("--out", required=True, help=".stat or .stat.gz")
    ap.add_argument("--summary", help="write summary JSON here")
    ap.add_argument("--chunk", type=int, default=5_000_000)
    args = ap.parse_args()

    n_records = 0
    parts = []
    reader = pd.read_csv(sys.stdin.buffer, sep=" ", header=None, names=["a", "b", "c", "bin"],
                         dtype=np.int64, chunksize=args.chunk, engine="c")
    for chunk in reader:
        n_records += len(chunk)
        if (chunk[["a", "b", "c"]].to_numpy() >= 1 << 20).any() or (chunk.to_numpy() < 0).any():
            sys.exit("barcode component out of range")
        code = (chunk["a"].to_numpy() << 40) | (chunk["b"].to_numpy() << 20) | chunk["c"].to_numpy()
        pairs = np.unique(np.stack([code, chunk["bin"].to_numpy()], axis=1), axis=0)
        parts.append(pairs)
        if len(parts) >= 8:
            parts = [np.unique(np.concatenate(parts), axis=0)]
    if not parts:
        sys.exit("no records on stdin")
    pairs = np.unique(np.concatenate(parts), axis=0)

    codes, rows = np.unique(pairs[:, 0], return_inverse=True)
    bins = pairs[:, 1]
    n_bins = int(bins.max()) + 1
    B = sparse.csr_matrix((np.ones(len(bins), dtype=np.int32), (rows, bins)),
                          shape=(len(codes), n_bins))
    M = (B.T @ B).tocoo()
    order = np.lexsort((M.col, M.row))
    i, j, c = M.row[order].astype(np.int64), M.col[order].astype(np.int64), M.data[order].astype(np.int64)

    opener = gzip.open if args.out.endswith(".gz") else open
    with opener(args.out, "wt", compresslevel=1) if args.out.endswith(".gz") else open(args.out, "w") as fh:
        pd.DataFrame({"a": i * args.bin_size, "b": j * args.bin_size, "c": c}).to_csv(
            fh, sep="\t", header=False, index=False)

    diag = i == j
    summary = {
        "records": int(n_records), "barcode_bin_pairs": int(len(bins)), "barcodes": int(len(codes)),
        "n_bins": n_bins, "lines": int(len(c)), "upper_nnz": int((i <= j).sum()),
        "total": int(c.sum()), "diag": int(c[diag].sum()), "max": int(c.max()),
    }
    text = json.dumps(summary)
    print(text, file=sys.stderr)
    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()

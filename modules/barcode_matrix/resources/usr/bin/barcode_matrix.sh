#!/usr/bin/env bash
# Barcode co-occurrence table for one chromosome of a duplicate-marked BAM whose read
# names end in "#A_B_C" (the stLFR barcode).
#
#   barcode_matrix.sh BAM CHROM BIN_SIZE OUT.stat.gz SUMMARY.json [SAMTOOLS_THREADS]
#
# Reads with MAPQ < 20, duplicates (flag 1024) and the null barcode 0_0_0 are dropped.
# Each remaining alignment record contributes (barcode, int(POS / BIN_SIZE)).
set -euo pipefail
BAM=$1
CHROM=$2
BIN=$3
OUT=$4
SUMMARY=$5
THREADS=${6:-2}

samtools view -@ "$THREADS" -F 1024 -q 20 "$BAM" "$CHROM" \
  | awk -F'\t' -v bin=$BIN '{n=split($1,a,"#"); if(n<2){bad++; next} bc=a[n]; if(bc=="0_0_0") next; gsub("_"," ",bc); print bc, int($4/bin)} END{if(bad) print bad" records without barcode" > "/dev/stderr"}' \
  | barcode_share_matrix.py --bin-size "$BIN" --out "$OUT" --summary "$SUMMARY"

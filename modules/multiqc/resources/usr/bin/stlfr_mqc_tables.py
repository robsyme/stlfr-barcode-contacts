#!/usr/bin/env python3
"""Turn the pipeline's per-sample statistics into MultiQC custom-content tables.

  stlfr_mqc_tables.py DIR

Reads from DIR and writes next to them:
  <sample>.tag_stats.json                    -> barcode_tagging_mqc.tsv
  <sample>.<chrom>_<bin>_bin.summary.json    -> barcode_matrix_mqc.tsv

A table is skipped when none of its input files is present. Standard library only.
samtools markdup, flagstat and stats files need no conversion: MultiQC reads them itself.
"""
import json
import re
import sys
from pathlib import Path


def write_table(path, header, columns, rows):
    """columns: list of (key, title, description, format). rows: list of (name, dict)."""
    lines = [f"# {k}: {json.dumps(v)}" for k, v in header.items()]
    lines.append("# headers:")
    for key, title, desc, fmt in columns:
        lines.append(f"#     {key}:")
        lines.append(f"#         title: {json.dumps(title)}")
        lines.append(f"#         description: {json.dumps(desc)}")
        lines.append(f"#         format: {json.dumps(fmt)}")
    lines.append("\t".join(["Sample"] + [c[0] for c in columns]))
    for name, values in rows:
        lines.append("\t".join([name] + [str(values.get(c[0], "")) for c in columns]))
    Path(path).write_text("\n".join(lines) + "\n")


def tagging(d):
    rows = []
    for f in sorted(d.glob("*.tag_stats.json")):
        s = json.loads(f.read_text())
        rows.append((f.name[: -len(".tag_stats.json")], {
            "total_pairs": s.get("total_pairs", ""),
            "valid_pct": round(100 * s["valid_fraction"], 2) if "valid_fraction" in s else "",
            "kept_pairs": s.get("kept_pairs", ""),
            "kept_pct": round(100 * s["kept_pairs"] / s["total_pairs"], 2) if s.get("total_pairs") else "",
            "barcodes_seen": s.get("distinct_barcodes_seen", ""),
            "barcodes_kept": s.get("distinct_barcodes_kept", ""),
            "barcode_fraction": s.get("barcode_fraction", ""),
        }))
    if rows:
        write_table(d / "barcode_tagging_mqc.tsv", {
            "id": "barcode_tagging",
            "section_name": "Barcode tagging",
            "description": "stLFR barcodes read from the end of read 2. A pair is valid when all three barcodes "
                           "match the whitelist with at most one substitution. Kept pairs are the valid pairs "
                           "whose barcode survived the barcode subsample.",
            "plot_type": "table",
        }, [
            ("total_pairs", "Read pairs", "Read pairs in the input", "{:,.0f}"),
            ("valid_pct", "% valid", "Percentage of read pairs with a valid barcode", "{:,.2f}"),
            ("kept_pairs", "Kept pairs", "Read pairs written for alignment", "{:,.0f}"),
            ("kept_pct", "% kept", "Kept pairs as a percentage of all read pairs", "{:,.2f}"),
            ("barcodes_seen", "Barcodes seen", "Distinct valid barcodes in the input", "{:,.0f}"),
            ("barcodes_kept", "Barcodes kept", "Distinct barcodes after the barcode subsample", "{:,.0f}"),
            ("barcode_fraction", "Barcode fraction", "Requested fraction of barcodes to keep", "{:,.2f}"),
        ], rows)
    return len(rows)


def natural(text):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", text)]


def matrix(d):
    rows = []
    for f in sorted(d.glob("*_bin.summary.json"), key=lambda p: natural(p.name)):
        m = re.fullmatch(r"(.+)_([^_]+)_bin\.summary\.json", f.name)
        if not m:
            continue
        s = json.loads(f.read_text())
        barcodes = s.get("barcodes", 0)
        rows.append((m.group(1), {
            "bin": m.group(2),
            "records": s.get("records", ""),
            "barcodes": s.get("barcodes", ""),
            "bins_per_barcode": round(s["barcode_bin_pairs"] / barcodes, 3) if barcodes and "barcode_bin_pairs" in s else "",
            "bin_pairs": s.get("upper_nnz", ""),
            "total": s.get("total", ""),
            "diag": s.get("diag", ""),
            "max": s.get("max", ""),
        }))
    if rows:
        write_table(d / "barcode_matrix_mqc.tsv", {
            "id": "barcode_matrix",
            "section_name": "Barcode co-occurrence tables",
            "description": "One row per sample and chromosome. Counts are distinct barcodes shared by a pair "
                           "of bins; the diagonal holds the barcodes seen in each bin.",
            "plot_type": "table",
        }, [
            ("bin", "Bin size", "Bin size of the table", "{}"),
            ("records", "Alignments", "Alignment records used (MAPQ >= 20, not duplicate, barcode not 0_0_0)", "{:,.0f}"),
            ("barcodes", "Barcodes", "Distinct barcodes on the chromosome", "{:,.0f}"),
            ("bins_per_barcode", "Bins per barcode", "Mean number of bins a barcode is seen in", "{:,.3f}"),
            ("bin_pairs", "Bin pairs", "Bin pairs with a count, diagonal included, each pair once", "{:,.0f}"),
            ("total", "Sum of counts", "Sum over all lines of the table (both orientations)", "{:,.0f}"),
            ("diag", "Diagonal sum", "Sum of the diagonal counts", "{:,.0f}"),
            ("max", "Max count", "Largest count in the table", "{:,.0f}"),
        ], rows)
    return len(rows)


def main():
    d = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    print(f"tagging rows: {tagging(d)}, matrix rows: {matrix(d)}", file=sys.stderr)


if __name__ == "__main__":
    main()

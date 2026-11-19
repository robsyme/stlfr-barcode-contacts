# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///
"""Tests for modules/tag_barcodes/resources/usr/bin/stlfr_tag_barcodes.py.

Run with either of:
    uv run --with pytest pytest tests/test_stlfr_tag_barcodes.py
    uv run tests/test_stlfr_tag_barcodes.py

All read data is synthetic and generated here.
"""
import gzip
import importlib.util
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RESOURCES = REPO / "modules" / "tag_barcodes" / "resources"
SCRIPT = RESOURCES / "usr" / "bin" / "stlfr_tag_barcodes.py"
WHITELIST = RESOURCES / "barcode_whitelist.tsv"

_spec = importlib.util.spec_from_file_location("stlfr_tag_barcodes", SCRIPT)
tag = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tag)

ENTRIES = tag.load_whitelist(WHITELIST)
ID_OF = dict(ENTRIES)
SEQS = [seq for seq, _ in ENTRIES]


# --- synthetic data ----------------------------------------------------------

def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng, seq, n):
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), n):
        seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
    return "".join(seq)


def make_pairs(n, seed=0, insert=100, spacer=6, n_molecules=400, style="sra"):
    """Synthetic stLFR pairs. Returns a list of dicts with the expected outcome.

    About 80% exact barcodes, 10% with one substitution in one barcode (still
    valid), 10% broken (two substitutions in one barcode, or an N).
    """
    rng = random.Random(seed)
    molecules = [tuple(rng.choice(SEQS) for _ in range(3)) for _ in range(n_molecules)]
    pairs = []
    for i in range(1, n + 1):
        bcs = list(rng.choice(molecules))
        expected = "_".join(str(ID_OF[b]) for b in bcs)
        roll = rng.random()
        which = rng.randrange(3)
        if roll < 0.10:
            bcs[which] = mutate(rng, bcs[which], 1)
        elif roll < 0.17:
            bcs[which] = mutate(rng, bcs[which], 2)
            # Two substitutions can land on another valid barcode only by a
            # 1-in-many chance; recompute the truth instead of assuming.
            expected = None
        elif roll < 0.20:
            pos = rng.randrange(10)
            bcs[which] = bcs[which][:pos] + "N" + bcs[which][pos + 1:]
            expected = None
        s1 = rand_seq(rng, 100)
        ins = rand_seq(rng, insert)
        s2 = ins + bcs[0] + rand_seq(rng, spacer) + bcs[1] + rand_seq(rng, spacer) + bcs[2]
        q1 = "".join(rng.choice("FGHI:,") for _ in range(len(s1)))
        q2 = "".join(rng.choice("FGHI:,") for _ in range(len(s2)))
        if expected is None:
            ids = [BC_MAP.get(b.encode()) for b in bcs]
            expected = "_".join(map(str, ids)) if all(ids) else None
        if style == "sra":
            h1 = f"@SIM.{i} FLOWCELL{i:08d} length={len(s1)}"
            h2 = f"@SIM.{i} FLOWCELL{i:08d} length={len(s2)}"
            plus1, plus2 = "+" + h1[1:], "+" + h2[1:]
        else:
            h1, h2, plus1, plus2 = f"@SIM.{i}/1", f"@SIM.{i}/2", "+", "+"
        pairs.append(dict(name=f"SIM.{i}", expected=expected, s1=s1, q1=q1, s2=s2, q2=q2,
                          insert=ins, rec1=f"{h1}\n{s1}\n{plus1}\n{q1}\n",
                          rec2=f"{h2}\n{s2}\n{plus2}\n{q2}\n"))
    return pairs


BC_MAP = tag.build_barcode_map(ENTRIES)


def write_fastq(path, text, gz=True):
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)


def write_pairs(tmp, pairs, gz=True):
    ext = ".fastq.gz" if gz else ".fastq"
    r1, r2 = tmp / f"in_1{ext}", tmp / f"in_2{ext}"
    write_fastq(r1, "".join(p["rec1"] for p in pairs), gz)
    write_fastq(r2, "".join(p["rec2"] for p in pairs), gz)
    return r1, r2


def run_cli(tmp, r1, r2, *extra, prefix="out/s", whitelist=WHITELIST):
    prefix = tmp / prefix
    cmd = [sys.executable, str(SCRIPT), "--r1", str(r1), "--r2", str(r2),
           "--whitelist", str(whitelist), "--out-prefix", str(prefix),
           "--stats", str(prefix) + ".tag_stats.json", *map(str, extra)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc, prefix


def read_fastq(path):
    with gzip.open(path, "rt") as fh:
        lines = fh.read().split("\n")
    assert lines[-1] == ""
    lines = lines[:-1]
    assert len(lines) % 4 == 0
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def read_chunks(prefix):
    chunks = []
    i = 0
    while True:
        f1 = Path(f"{prefix}.chunk{i:04d}_1.fastq.gz")
        f2 = Path(f"{prefix}.chunk{i:04d}_2.fastq.gz")
        if not f1.exists():
            assert not f2.exists()
            break
        chunks.append((read_fastq(f1), read_fastq(f2)))
        i += 1
    return chunks


def stats_of(prefix):
    return json.loads(Path(str(prefix) + ".tag_stats.json").read_text())


# --- unit tests --------------------------------------------------------------

def test_whitelist_is_bgi_numbering():
    assert len(ENTRIES) == 1536
    assert [idx for _, idx in ENTRIES] == list(range(1, 1537))
    assert ENTRIES[0] == ("TAACAGCCAA", 1)


def test_barcode_map_exact_and_one_substitution():
    seq, idx = ENTRIES[41]
    assert BC_MAP[seq.encode()] == idx
    rng = random.Random(1)
    for _ in range(200):
        seq, idx = rng.choice(ENTRIES)
        assert BC_MAP[mutate(rng, seq, 1).encode()] == idx
    assert BC_MAP.get(("N" + seq[1:]).encode()) is None
    assert BC_MAP.get(seq[:9].encode()) is None
    assert len(BC_MAP) == 1536 * 31


def test_barcode_map_later_entry_wins_like_split_barcode_pl():
    # AAAAAAAAAA and AAAAAAAAAC share the variants AAAAAAAAA[ACGT].
    entries = [("AAAAAAAAAA", 1), ("AAAAAAAAAC", 2)]
    m = tag.build_barcode_map(entries)
    assert m[b"AAAAAAAAAA"] == 2
    assert m[b"AAAAAAAAAG"] == 2
    assert m[b"CAAAAAAAAA"] == 1


def test_whitelist_without_id_column(tmp_path):
    path = tmp_path / "wl.txt"
    path.write_text("".join(seq + "\n" for seq in SEQS[:50]))
    assert tag.load_whitelist(path) == ENTRIES[:50]


def test_subsample_hash_is_pinned():
    # Golden values: these must never change, on any platform or Python.
    tables = tag.hash_tables(1, 1537)
    assert tables[0][1] == 18215591310788271736
    assert tables[2][1536] == 8851242830403249417
    threshold = tag.keep_threshold(0.25)
    got = [tag.barcode_is_kept(a, b, c, tables, threshold) for a, b, c in
           [(1, 1, 1), (1283, 1147, 136), (710, 226, 170), (1536, 1536, 1536),
            (5, 6, 7), (7, 6, 5), (100, 200, 300), (2, 2, 2)]]
    assert got == [False, True, False, False, False, False, False, False]


@pytest.mark.parametrize("fraction", [0.05, 0.25, 0.5, 0.9])
def test_subsample_fraction_is_close(fraction):
    tables = tag.hash_tables(7, 1537)
    threshold = tag.keep_threshold(fraction)
    rng = random.Random(3)
    n = 40_000
    kept = sum(tag.barcode_is_kept(rng.randint(1, 1536), rng.randint(1, 1536),
                                   rng.randint(1, 1536), tables, threshold) for _ in range(n))
    assert abs(kept / n - fraction) < 0.01


def test_subsample_fraction_one_keeps_everything():
    tables = tag.hash_tables(1, 1537)
    threshold = tag.keep_threshold(1.0)
    assert all(tag.barcode_is_kept(a, a, a, tables, threshold) for a in range(1, 1537))


def test_cut_after_lines():
    rng = random.Random(5)
    for _ in range(200):
        lines = [b"x" * rng.choice([0, 1, 5, 150, 150, 150, 4000]) for _ in range(rng.randint(1, 300))]
        buf = b"".join(line + b"\n" for line in lines)
        for n in {0, 1, len(lines) // 2, len(lines) - 1, len(lines)}:
            cut = tag.cut_after_lines(buf, n, len(lines))
            assert buf[:cut] == b"".join(line + b"\n" for line in lines[:n])


def test_detect_offsets_for_other_layouts():
    for insert, spacer in [(100, 6), (90, 6), (100, 0), (126, 6)]:
        pairs = make_pairs(500, seed=insert + spacer, insert=insert, spacer=spacer)
        offsets, rates = tag.detect_barcode_offsets([p["s2"].encode() for p in pairs], BC_MAP)
        assert offsets == [insert, insert + 10 + spacer, insert + 20 + 2 * spacer]
        assert all(rate > 0.8 for rate in rates.values())


# --- end to end --------------------------------------------------------------

@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("data")
    pairs = make_pairs(3000, seed=11)
    r1, r2 = write_pairs(tmp, pairs)
    return pairs, r1, r2


def check_against_truth(pairs, chunks, keep=lambda expected: True):
    out1 = [rec for c in chunks for rec in c[0]]
    out2 = [rec for c in chunks for rec in c[1]]
    want = [p for p in pairs if p["expected"] and keep(p["expected"])]
    assert len(out1) == len(out2) == len(want)
    for p, rec1, rec2 in zip(want, out1, out2):
        assert rec1 == (f"@{p['name']}#{p['expected']}/1", p["s1"], "+", p["q1"])
        assert rec2 == (f"@{p['name']}#{p['expected']}/2", p["insert"], "+",
                        p["q2"][:len(p["insert"])])
    return want


def test_end_to_end_all_barcodes(tmp_path, dataset):
    pairs, r1, r2 = dataset
    proc, prefix = run_cli(tmp_path, r1, r2, "--threads", 1)
    assert proc.returncode == 0, proc.stderr
    chunks = read_chunks(prefix)
    assert len(chunks) == 1
    want = check_against_truth(pairs, chunks)
    n_valid = len(want)
    assert 0.8 * len(pairs) < n_valid < len(pairs)
    stats = stats_of(prefix)
    assert stats["total_pairs"] == 3000
    assert stats["valid_barcode_pairs"] == stats["kept_pairs"] == n_valid
    assert stats["dropped_invalid"] == 3000 - n_valid
    assert stats["dropped_by_subsample"] == 0
    n_distinct = len({p["expected"] for p in want})
    assert stats["distinct_barcodes_seen"] == stats["distinct_barcodes_kept"] == n_distinct
    assert stats["barcode_offsets"] == [100, 116, 132]
    assert stats["n_chunks"] == 1
    assert stats["chunk_pairs"] == [n_valid]
    assert stats["chunk_files"] == [["s.chunk0000_1.fastq.gz", "s.chunk0000_2.fastq.gz"]]
    assert (stats["r1_len_mode"], stats["r2_len_mode"]) == (100, 142)


def test_chunking_and_threads_give_identical_reads(tmp_path, dataset):
    pairs, r1, r2 = dataset
    results = []
    for name, extra in [("a", ["--threads", 1]),
                        ("b", ["--threads", 4, "--block-mb", 0.05]),
                        ("c", ["--threads", 3, "--block-mb", 0.011])]:
        proc, prefix = run_cli(tmp_path, r1, r2, "--chunk-pairs", 700, *extra,
                               prefix=f"{name}/s")
        assert proc.returncode == 0, proc.stderr
        chunks = read_chunks(prefix)
        n_valid = len(check_against_truth(pairs, chunks))
        sizes = [len(c[0]) for c in chunks]
        assert sizes == [700] * (n_valid // 700) + ([n_valid % 700] if n_valid % 700 else [])
        stats = stats_of(prefix)
        assert stats["n_chunks"] == len(chunks)
        assert stats["chunk_pairs"] == sizes
        results.append((chunks, {k: stats[k] for k in (
            "total_pairs", "valid_barcode_pairs", "kept_pairs", "distinct_barcodes_seen",
            "distinct_barcodes_kept")}))
    assert results[0] == results[1] == results[2]


def test_exact_multiple_of_chunk_size_leaves_no_empty_chunk(tmp_path, dataset):
    pairs, r1, r2 = dataset
    n_valid = sum(1 for p in pairs if p["expected"])
    proc, prefix = run_cli(tmp_path, r1, r2, "--chunk-pairs", n_valid, "--threads", 2,
                           "--block-mb", 0.05)
    assert proc.returncode == 0, proc.stderr
    assert [len(c[0]) for c in read_chunks(prefix)] == [n_valid]


def test_barcode_fraction_keeps_whole_barcodes(tmp_path, dataset):
    pairs, r1, r2 = dataset
    tables = tag.hash_tables(1, 1537)
    threshold = tag.keep_threshold(0.25)

    def keep(expected):
        return tag.barcode_is_kept(*map(int, expected.split("_")), tables, threshold)

    outputs = []
    for name, threads in [("a", 1), ("b", 4)]:
        proc, prefix = run_cli(tmp_path, r1, r2, "--barcode-fraction", 0.25, "--seed", 1,
                               "--threads", threads, "--block-mb", 0.05, prefix=f"{name}/s")
        assert proc.returncode == 0, proc.stderr
        chunks = read_chunks(prefix)
        want = check_against_truth(pairs, chunks, keep)
        stats = stats_of(prefix)
        valid = [p for p in pairs if p["expected"]]
        assert stats["kept_pairs"] == len(want)
        assert stats["dropped_by_subsample"] == len(valid) - len(want)
        assert stats["distinct_barcodes_seen"] == len({p["expected"] for p in valid})
        assert stats["distinct_barcodes_kept"] == len({p["expected"] for p in want})
        assert 0.15 < stats["distinct_barcodes_kept"] / stats["distinct_barcodes_seen"] < 0.35
        outputs.append(chunks)
    assert outputs[0] == outputs[1]

    proc, prefix = run_cli(tmp_path, r1, r2, "--barcode-fraction", 0.25, "--seed", 2,
                           "--threads", 1, prefix="seed2/s")
    assert proc.returncode == 0, proc.stderr
    assert read_chunks(prefix) != outputs[0]


def test_mate_suffixes_and_plain_fastq(tmp_path):
    pairs = make_pairs(300, seed=4, style="slash")
    r1, r2 = write_pairs(tmp_path, pairs, gz=False)
    proc, prefix = run_cli(tmp_path, r1, r2, "--threads", 1)
    assert proc.returncode == 0, proc.stderr
    check_against_truth(pairs, read_chunks(prefix))


def test_other_insert_length_is_trimmed_at_detected_offset(tmp_path):
    pairs = make_pairs(400, seed=9, insert=90)
    r1, r2 = write_pairs(tmp_path, pairs)
    proc, prefix = run_cli(tmp_path, r1, r2, "--threads", 1)
    assert proc.returncode == 0, proc.stderr
    check_against_truth(pairs, read_chunks(prefix))
    assert stats_of(prefix)["barcode_offsets"] == [90, 106, 122]


def test_final_newline_missing_is_accepted(tmp_path):
    pairs = make_pairs(200, seed=2)
    write_fastq(tmp_path / "a_1.fastq", "".join(p["rec1"] for p in pairs)[:-1], gz=False)
    write_fastq(tmp_path / "a_2.fastq", "".join(p["rec2"] for p in pairs)[:-1], gz=False)
    proc, prefix = run_cli(tmp_path, tmp_path / "a_1.fastq", tmp_path / "a_2.fastq",
                           "--threads", 1)
    assert proc.returncode == 0, proc.stderr
    check_against_truth(pairs, read_chunks(prefix))


# --- failures ----------------------------------------------------------------

def fails_with(proc, *needles):
    assert proc.returncode != 0
    assert "Traceback" not in proc.stderr, proc.stderr
    for needle in needles:
        assert needle in proc.stderr, proc.stderr


@pytest.mark.parametrize("threads", [1, 3])
def test_r2_has_fewer_pairs(tmp_path, threads):
    pairs = make_pairs(600, seed=1)
    r1, _ = write_pairs(tmp_path, pairs)
    r2 = tmp_path / "short_2.fastq.gz"
    write_fastq(r2, "".join(p["rec2"] for p in pairs[:450]))
    proc, _ = run_cli(tmp_path, r1, r2, "--threads", threads, "--block-mb", 0.02)
    fails_with(proc, "different numbers of reads", "R2 ends after 450")


def test_r1_has_fewer_pairs(tmp_path):
    pairs = make_pairs(600, seed=1)
    _, r2 = write_pairs(tmp_path, pairs)
    r1 = tmp_path / "short_1.fastq.gz"
    write_fastq(r1, "".join(p["rec1"] for p in pairs[:450]))
    proc, _ = run_cli(tmp_path, r1, r2, "--threads", 2, "--block-mb", 0.02)
    fails_with(proc, "different numbers of reads", "R1 ends after 450")


@pytest.mark.parametrize("mate", [1, 2])
def test_truncated_mid_record(tmp_path, mate):
    pairs = make_pairs(300, seed=1)
    texts = ["".join(p[f"rec{m}"] for p in pairs) for m in (1, 2)]
    texts[mate - 1] = texts[mate - 1][:len(texts[mate - 1]) - 150]
    write_fastq(tmp_path / "t_1.fastq.gz", texts[0])
    write_fastq(tmp_path / "t_2.fastq.gz", texts[1])
    proc, _ = run_cli(tmp_path, tmp_path / "t_1.fastq.gz", tmp_path / "t_2.fastq.gz",
                      "--threads", 1)
    fails_with(proc, "truncated")


def test_truncated_gzip_stream(tmp_path):
    pairs = make_pairs(2000, seed=1)
    r1, r2 = write_pairs(tmp_path, pairs)
    data = r2.read_bytes()
    cut = tmp_path / "cut_2.fastq.gz"
    cut.write_bytes(data[:len(data) * 2 // 3])
    proc, _ = run_cli(tmp_path, r1, cut, "--threads", 2)
    fails_with(proc)
    assert "truncated" in proc.stderr.lower() or "different numbers" in proc.stderr


def test_out_of_sync_names(tmp_path):
    pairs = make_pairs(400, seed=1)
    r1, _ = write_pairs(tmp_path, pairs)
    shuffled = pairs[:200] + [pairs[201], pairs[200]] + pairs[202:]
    r2 = tmp_path / "swap_2.fastq.gz"
    write_fastq(r2, "".join(p["rec2"] for p in shuffled))
    proc, _ = run_cli(tmp_path, r1, r2, "--threads", 2)
    fails_with(proc, "out of sync", "pair 201", "SIM.201", "SIM.202")


def test_name_prefix_is_not_mistaken_for_a_match(tmp_path):
    # SIM.1 is a prefix of SIM.10: the check must compare whole names.
    pairs = make_pairs(20, seed=1)
    r1, _ = write_pairs(tmp_path, pairs)
    wrong = [pairs[9]] + pairs[1:]
    r2 = tmp_path / "prefix_2.fastq.gz"
    write_fastq(r2, "".join(p["rec2"] for p in wrong))
    proc, _ = run_cli(tmp_path, r1, r2, "--threads", 1)
    fails_with(proc, "out of sync", "pair 1:")


def test_low_valid_fraction_fails(tmp_path):
    rng = random.Random(8)
    pairs = make_pairs(500, seed=1)
    text2 = ""
    for i, p in enumerate(pairs):
        s2 = p["s2"] if i % 4 == 0 else p["insert"] + rand_seq(rng, 42)
        text2 += f"@{p['name']}\n{s2}\n+\n{p['q2']}\n"
    r1, _ = write_pairs(tmp_path, pairs)
    r2 = tmp_path / "low_2.fastq.gz"
    write_fastq(r2, text2)
    proc, _ = run_cli(tmp_path, r1, r2, "--threads", 1)
    fails_with(proc, "valid", "minimum 50%")
    proc, prefix = run_cli(tmp_path, r1, r2, "--threads", 1, "--min-valid-fraction", 0.1,
                           prefix="lenient/s")
    assert proc.returncode == 0, proc.stderr
    assert 0.15 < stats_of(prefix)["valid_fraction"] < 0.35


def test_bad_fraction_is_rejected(tmp_path, dataset):
    _, r1, r2 = dataset
    for value in (0, 1.5):
        proc, _ = run_cli(tmp_path, r1, r2, "--barcode-fraction", value)
        fails_with(proc, "--barcode-fraction")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))

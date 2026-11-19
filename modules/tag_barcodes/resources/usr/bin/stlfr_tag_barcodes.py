#!/usr/bin/env python3
"""Assign stLFR barcodes to read pairs and write FASTQ chunks for alignment.

stLFR (MGI co-barcoded) libraries carry three 10 bp barcodes at the end of
read 2, after the genomic insert. This script does the job of BGI's
`split_barcode.pl`: it looks each barcode up in the whitelist with one
substitution tolerated, renames both reads to `name#A_B_C/1` and `/2`, trims
read 2 back to the insert, and drops pairs without a valid barcode (BGI's
`0_0_0`). It can also keep a deterministic subset of barcodes.

Barcode numbering and matching follow `split_barcode.pl`: A, B and C are the
IDs from column 2 of `barcode.list`, and when a one-substitution variant could
belong to two whitelist entries the later entry in the file wins.

The barcode offsets are detected from the data, since deposited read 2 lengths
do not always match the textbook layouts.

Standard library only. Decompression runs in external `pigz` (or `gzip`)
processes; parsing and gzip compression run in `--threads` worker processes.
Output files are multi-member gzip, which every gzip reader accepts.
"""
import argparse
import gzip
import json
import multiprocessing
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
import zlib
from array import array
from collections import Counter, deque

BC_LEN = 10
N_BARCODES = 3
_BASES = "ACGT"
_MASK64 = (1 << 64) - 1
_NAME_END = (b"", b" ", b"\t", b"/")
_MATE_SUFFIX = (b"/1", b"/2")


class InputError(Exception):
    """Malformed, truncated or out-of-sync input, or a failed quality gate."""


# --- whitelist and matching --------------------------------------------------

def load_whitelist(path):
    """Read `SEQ [ID]` lines. Returns a list of (sequence, integer ID).

    Without an ID column the ID is the 1-based line number, which is how BGI's
    `barcode.list` is numbered anyway.
    """
    entries = []
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            seq = fields[0].upper()
            if len(seq) != BC_LEN or set(seq) - set(_BASES):
                raise InputError(f"whitelist {path}: bad barcode sequence {fields[0]!r}")
            idx = int(fields[1]) if len(fields) > 1 else len(entries) + 1
            if idx < 1:
                raise InputError(f"whitelist {path}: barcode IDs must be >= 1, got {idx}")
            entries.append((seq, idx))
    if not entries:
        raise InputError(f"whitelist {path} is empty")
    return entries


def build_barcode_map(entries):
    """Map every whitelist sequence and one-substitution variant to its ID.

    Built the way `split_barcode.pl` builds its hash: entries in file order,
    every position replaced by each of A/C/G/T, later assignments overwriting
    earlier ones. Keys are bytes.
    """
    bc_map = {}
    for seq, idx in entries:
        for pos in range(len(seq)):
            for base in _BASES:
                bc_map[(seq[:pos] + base + seq[pos + 1:]).encode()] = idx
    return bc_map


# --- barcode subsampling -----------------------------------------------------

def _splitmix64(state):
    state = (state + 0x9E3779B97F4A7C15) & _MASK64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
    return state, z ^ (z >> 31)


def hash_tables(seed, radix):
    """Three tables of 64-bit values, one per barcode position.

    The hash of a barcode is `t1[A] ^ t2[B] ^ t3[C]` (tabulation hashing). The
    tables come from a splitmix64 stream seeded with `seed`, so the result
    depends only on (A, B, C, seed), on any machine and Python version.
    """
    state = (int(seed) * 0x9E3779B97F4A7C15 + 0x632BE59BD9B4E019) & _MASK64
    tables = []
    for _ in range(N_BARCODES):
        table = []
        for _ in range(radix):
            state, value = _splitmix64(state)
            table.append(value)
        tables.append(table)
    return tables


def keep_threshold(fraction):
    return int(fraction * (1 << 64))


def barcode_is_kept(a, b, c, tables, threshold):
    return (tables[0][a] ^ tables[1][b] ^ tables[2][c]) < threshold


# --- offset detection --------------------------------------------------------

def detect_barcode_offsets(r2_seqs, bc_map, n_barcodes=N_BARCODES, bc_len=BC_LEN):
    """Find the barcode start positions in read 2.

    Scores every window of the modal-length reads by its whitelist match rate
    and greedily takes the best non-overlapping windows. Returns
    (sorted offsets, {offset: match rate}).
    """
    if not r2_seqs:
        raise InputError("no reads available to detect barcode offsets")
    mode_len = Counter(len(s) for s in r2_seqs).most_common(1)[0][0]
    seqs = [s for s in r2_seqs if len(s) == mode_len]
    total = len(seqs)
    rates = []
    for start in range(0, mode_len - bc_len + 1):
        end = start + bc_len
        hits = sum(1 for s in seqs if s[start:end] in bc_map)
        rates.append((hits / total, start))
    rates.sort(key=lambda x: (-x[0], x[1]))
    chosen = {}
    for rate, start in rates:
        if all(abs(start - c) >= bc_len for c in chosen):
            chosen[start] = rate
        if len(chosen) == n_barcodes:
            break
    if len(chosen) < n_barcodes:
        raise InputError(
            f"read 2 (length {mode_len}) is too short to hold {n_barcodes} barcodes")
    return sorted(chosen), {k: chosen[k] for k in sorted(chosen)}


# --- input -------------------------------------------------------------------

def _decompressor():
    for tool in ("pigz", "gzip"):
        path = shutil.which(tool)
        if path:
            return [path, "-dc"]
    return None


class Source:
    """A FASTQ file read in the background, decompressed externally.

    A thread pulls pieces off the decompressor into a bounded queue, so that
    the R1 and R2 decompressors keep running while the main loop is busy with
    the other file.
    """

    def __init__(self, path, piece=4 << 20):
        self.path = path
        self.piece = piece
        self.proc = None
        with open(path, "rb") as fh:
            gzipped = fh.read(2) == b"\x1f\x8b"
        cmd = _decompressor() if gzipped else None
        if not gzipped:
            self.fh = open(path, "rb", buffering=1 << 20)
        elif cmd is None:
            self.fh = gzip.open(path, "rb")
        else:
            self.proc = subprocess.Popen(cmd + [path], stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, bufsize=1 << 20)
            self.fh = self.proc.stdout
        self.queue = queue.Queue(maxsize=8)
        self.done = False
        self.thread = threading.Thread(target=self._pump, daemon=True)
        self.thread.start()

    def _pump(self):
        try:
            while True:
                data = self.fh.read(self.piece)
                self.queue.put(data)
                if not data:
                    return
        except Exception as exc:  # handed to the main thread by read()
            self.queue.put(exc)

    def read(self):
        """The next piece of the file; b"" at the end, and from then on."""
        if self.done:
            return b""
        data = self.queue.get()
        if isinstance(data, Exception):
            self.done = True
            raise InputError(f"{self.path}: cannot decompress ({data})")
        if not data:
            self.done = True
        return data

    def finish(self):
        """Call at end of input. Raises if the decompressor reported an error."""
        if self.proc is None:
            self.fh.close()
            return
        err = self.proc.stderr.read().decode(errors="replace").strip()
        rc = self.proc.wait()
        if rc != 0:
            raise InputError(
                f"{self.path}: decompression failed (exit {rc}): {err or 'no message'}; "
                "the file is probably truncated or corrupt")

    def kill(self):
        if self.proc is not None and self.proc.poll() is None:
            self.proc.kill()
            self.proc.wait()


def cut_after_lines(buf, n, total_lines):
    """Offset just past the n-th newline of buf, which holds total_lines >= n.

    Jumps to the proportional position, counts newlines up to it (a fast C
    scan) and walks the last few lines. FASTQ lines are regular enough that the
    walk is short.
    """
    if n <= 0:
        return 0
    size = len(buf)
    pos = min(size, size * n // total_lines)
    count = buf.count(b"\n", 0, pos)
    for _ in range(8):
        diff = n - count
        if -64 <= diff <= 64:
            break
        new_pos = max(0, min(size, pos + diff * size // total_lines))
        if new_pos > pos:
            count += buf.count(b"\n", pos, new_pos)
        else:
            count -= buf.count(b"\n", new_pos, pos)
        pos = new_pos
    diff = n - count
    if diff > 0:
        for _ in range(diff):
            pos = buf.index(b"\n", pos) + 1
        return pos
    for _ in range(1 - diff):
        pos = buf.rfind(b"\n", 0, pos)
    return pos + 1


class BlockReader:
    """Hands out blocks of whole lines from a Source."""

    def __init__(self, source, label):
        self.source = source
        self.label = label
        self.buf = b""
        self.lines = 0
        self.eof = False

    def _fill(self):
        data = self.source.read()
        if not data:
            self.eof = True
            if self.buf and not self.buf.endswith(b"\n"):
                self.buf += b"\n"
                self.lines += 1
            return
        self.lines += data.count(b"\n")
        self.buf = self.buf + data if self.buf else data

    def take_records(self, target_bytes):
        """About target_bytes of whole 4-line records. Returns (bytes, n_lines);
        (b"", 0) at end of input."""
        while not self.eof and (len(self.buf) < target_bytes or self.lines < 4):
            self._fill()
        if self.eof and self.lines % 4:
            raise InputError(
                f"{self.label} is truncated: it ends part-way through a FASTQ record "
                f"({self.lines % 4} of 4 lines)")
        n = self.lines - self.lines % 4
        if len(self.buf) > 1.5 * target_bytes:
            # More is buffered than asked for: hand out a proportional share.
            n = max(4, int(self.lines * target_bytes / len(self.buf)) // 4 * 4)
        return self.take_lines(n), n

    def take_lines(self, n):
        """Exactly n lines, or fewer if the input ends first."""
        if n <= 0:
            return b""
        while self.lines < n and not self.eof:
            self._fill()
        if self.lines < n or (self.lines == n and self.buf.endswith(b"\n")):
            out, self.buf, self.lines = self.buf, b"", 0
            return out
        if self.lines == n:
            # The buffer ends with the start of the next line; leave it behind.
            cut = self.buf.rfind(b"\n") + 1
        else:
            cut = cut_after_lines(self.buf, n, self.lines)
        out, self.buf = self.buf[:cut], self.buf[cut:]
        self.lines -= n
        return out

    def at_end(self):
        if not self.buf and not self.eof:
            self._fill()
        return self.eof and not self.buf


# --- per-block work (runs in worker processes) --------------------------------

_W = {}


def init_worker(whitelist, fraction, seed, level):
    entries = load_whitelist(whitelist)
    radix = max(idx for _seq, idx in entries) + 1
    _W["bc_map"] = build_barcode_map(entries)
    _W["radix"] = radix
    _W["names"] = [b"%d" % i for i in range(radix)]
    _W["level"] = level
    if fraction < 1.0:
        _W["tables"] = hash_tables(seed, radix)
        _W["threshold"] = keep_threshold(fraction)
    else:
        _W["tables"] = None
        _W["threshold"] = None


def gzip_member(data, level):
    if not data:
        return b""
    comp = zlib.compressobj(level, zlib.DEFLATED, 31)
    return comp.compress(data) + comp.flush()


def process_block(job):
    """Tag one block of pairs.

    job is (index of the block's first pair, barcode offsets, R1 bytes, R2
    bytes); the two byte strings hold the same number of whole records. Returns (pairs, valid, kept, gzip R1, gzip
    R2, barcode keys as array('I') bytes, number of leading keys that are
    kept). A key is (A * radix + B) * radix + C; each distinct barcode of the
    block appears once, kept ones first.
    """
    first_pair, (o1, o2, o3), raw1, raw2 = job
    get = _W["bc_map"].get
    radix = _W["radix"]
    names = _W["names"]
    e1, e2, e3 = o1 + BC_LEN, o2 + BC_LEN, o3 + BC_LEN
    tables = _W["tables"]
    threshold = _W["threshold"]
    if tables is not None:
        t1, t2, t3 = tables

    l1 = raw1.split(b"\n")
    l2 = raw2.split(b"\n")
    del raw1, raw2
    n_lines = len(l1) - 1
    if len(l2) != len(l1):
        raise InputError("internal error: R1 and R2 blocks differ in size")

    out1 = []
    out2 = []
    add1 = out1.append
    add2 = out2.append
    kept_keys = set()
    dropped_keys = set()
    n_valid = 0
    n_kept = 0

    for i in range(0, n_lines, 4):
        h1 = l1[i]
        s2 = l2[i + 1]
        name = h1.split(None, 1)[0] if h1 else h1
        if name[-2:] in _MATE_SUFFIX:
            name = name[:-2]
        h2 = l2[i]
        s1 = l1[i + 1]
        q1 = l1[i + 3]
        q2 = l2[i + 3]
        # Structure and sync checks, on every pair.
        if (name[:1] != b"@" or not h2.startswith(name)
                or h2[len(name):len(name) + 1] not in _NAME_END
                or len(s1) != len(q1) or len(s2) != len(q2)
                or l1[i + 2][:1] != b"+" or l2[i + 2][:1] != b"+"):
            _explain_bad_pair(first_pair + i // 4, l1[i:i + 4], l2[i:i + 4])

        a = get(s2[o1:e1])
        if a is None:
            continue
        b = get(s2[o2:e2])
        if b is None:
            continue
        c = get(s2[o3:e3])
        if c is None:
            continue
        n_valid += 1
        key = (a * radix + b) * radix + c
        if tables is not None and (t1[a] ^ t2[b] ^ t3[c]) >= threshold:
            dropped_keys.add(key)
            continue
        kept_keys.add(key)
        n_kept += 1
        tag = b"#%b_%b_%b/" % (names[a], names[b], names[c])
        add1(b"%b%b1\n%b\n+\n%b\n" % (name, tag, s1, q1))
        add2(b"%b%b2\n%b\n+\n%b\n" % (name, tag, s2[:o1], q2[:o1]))

    level = _W["level"]
    gz1 = gzip_member(b"".join(out1), level)
    gz2 = gzip_member(b"".join(out2), level)
    keys = array("I", kept_keys)
    keys.extend(dropped_keys)
    return (n_lines // 4, n_valid, n_kept, gz1, gz2, keys.tobytes(), len(kept_keys))


def _explain_bad_pair(pair_index, rec1, rec2):
    where = f"pair {pair_index + 1}"

    def show(line):
        return line[:80].decode(errors="replace")

    h1, h2 = rec1[0], rec2[0]
    if h1[:1] != b"@" or h2[:1] != b"@":
        raise InputError(
            f"{where}: expected a FASTQ header starting with '@', found "
            f"R1 {show(h1)!r} and R2 {show(h2)!r}; the input is malformed or truncated")
    if rec1[2][:1] != b"+" or rec2[2][:1] != b"+":
        raise InputError(f"{where} ({show(h1)}): third line of the record is not '+'; "
                         "the input is malformed or truncated")
    if len(rec1[1]) != len(rec1[3]) or len(rec2[1]) != len(rec2[3]):
        raise InputError(f"{where} ({show(h1)}): sequence and quality lengths differ; "
                         "the input is malformed or truncated")
    raise InputError(
        f"{where}: R1 and R2 are out of sync: R1 read is {show(h1)!r} but R2 read is "
        f"{show(h2)!r}")


# --- output ------------------------------------------------------------------

class ChunkWriter:
    """Writes kept pairs to PREFIX.chunkNNNN_{1,2}.fastq.gz, chunk_pairs each."""

    def __init__(self, prefix, chunk_pairs, level):
        self.prefix = prefix
        self.chunk_pairs = chunk_pairs
        self.level = level
        self.n_chunks = 0
        self.in_chunk = 0
        self.fh = None
        self.files = []
        self.chunk_counts = []

    def _open(self):
        paths = [f"{self.prefix}.chunk{self.n_chunks:04d}_{mate}.fastq.gz" for mate in (1, 2)]
        self.fh = [open(p, "wb", buffering=1 << 20) for p in paths]
        self.files.append(paths)
        self.chunk_counts.append(0)
        self.n_chunks += 1
        self.in_chunk = 0

    def _close(self):
        if self.fh:
            for fh in self.fh:
                fh.close()
            self.fh = None

    def _write(self, gz1, gz2, n):
        if self.fh is None:
            self._open()
        self.fh[0].write(gz1)
        self.fh[1].write(gz2)
        self.in_chunk += n
        self.chunk_counts[-1] += n
        if self.in_chunk >= self.chunk_pairs:
            self._close()

    def write(self, gz1, gz2, n):
        """Append n pairs, given as one gzip member per mate."""
        if n == 0:
            return
        room = self.chunk_pairs - self.in_chunk if self.fh else self.chunk_pairs
        if n <= room:
            self._write(gz1, gz2, n)
            return
        # The block straddles a chunk boundary: split it on record boundaries.
        raw1 = zlib.decompress(gz1, 31)
        raw2 = zlib.decompress(gz2, 31)
        del gz1, gz2
        left = n
        while left:
            take = min(room, left)
            cut1 = cut_after_lines(raw1, 4 * take, 4 * left)
            cut2 = cut_after_lines(raw2, 4 * take, 4 * left)
            self._write(gzip_member(raw1[:cut1], self.level),
                        gzip_member(raw2[:cut2], self.level), take)
            raw1, raw2 = raw1[cut1:], raw2[cut2:]
            left -= take
            room = self.chunk_pairs

    def close(self):
        self._close()


class DistinctCounter:
    """Exact count of distinct barcodes, with one bit per possible barcode."""

    def __init__(self, radix):
        self.bits = bytearray(radix ** 3 // 8 + 1)
        self.seen = 0
        self.kept = 0

    def add(self, key_bytes, n_kept_keys):
        keys = array("I")
        keys.frombytes(key_bytes)
        bits = self.bits
        new = 0
        new_kept = 0
        for pos, key in enumerate(keys):
            byte = bits[key >> 3]
            bit = 1 << (key & 7)
            if not byte & bit:
                bits[key >> 3] = byte | bit
                new += 1
                if pos < n_kept_keys:
                    new_kept += 1
        self.seen += new
        self.kept += new_kept


# --- driver ------------------------------------------------------------------

def _first_seqs(block, limit):
    lines = block.split(b"\n", 4 * limit)
    return lines[1:4 * limit:4]


def run(args):
    t_start = time.time()
    if not 0 < args.barcode_fraction <= 1:
        raise InputError("--barcode-fraction must be in (0, 1]")
    if args.chunk_pairs < 1:
        raise InputError("--chunk-pairs must be at least 1")
    if not 0 <= args.compress_level <= 9:
        raise InputError("--compress-level must be 0-9")
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    entries = load_whitelist(args.whitelist)
    bc_map = build_barcode_map(entries)
    radix = max(idx for _seq, idx in entries) + 1
    if radix ** 3 > 1 << 32:
        raise InputError("whitelist IDs are too large (maximum ID is 1624)")

    threads = args.threads or len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") \
        else args.threads or os.cpu_count() or 1
    n_workers = max(1, threads - 2)

    # The pool is started before the reader threads, so that workers are never
    # forked from a multi-threaded process.
    init = (args.whitelist, args.barcode_fraction, args.seed, args.compress_level)
    block_bytes = max(1, int(args.block_mb * (1 << 20)))
    piece = max(1024, min(4 << 20, block_bytes // 2))
    pool = None
    sources = []
    writer = ChunkWriter(args.out_prefix, args.chunk_pairs, args.compress_level)
    try:
        if n_workers > 1:
            pool = multiprocessing.Pool(n_workers, initializer=init_worker, initargs=init)
            submit = lambda job: pool.apply_async(process_block, (job,))
            collect = lambda handle: handle.get()
        else:
            init_worker(*init)
            submit = process_block
            collect = lambda result: result
        sources.extend(Source(path, piece) for path in (args.r1, args.r2))
        r1 = BlockReader(sources[0], "R1")
        r2 = BlockReader(sources[1], "R2")

        def next_block():
            raw1, n_lines = r1.take_records(block_bytes)
            raw2 = r2.take_lines(n_lines)
            got = raw2.count(b"\n")
            if got != n_lines:
                raise InputError(
                    f"R1 and R2 hold different numbers of reads: R2 ends after "
                    f"{pairs_read + got // 4} complete pairs but R1 continues "
                    f"(R2 truncated?)" if got % 4 == 0 else
                    f"R2 is truncated: it ends part-way through a FASTQ record, "
                    f"after {pairs_read + got // 4} complete pairs")
            if n_lines == 0 and not r2.at_end():
                raise InputError(
                    f"R1 and R2 hold different numbers of reads: R1 ends after "
                    f"{pairs_read} pairs but R2 continues (R1 truncated?)")
            return raw1, raw2, n_lines // 4

        pairs_read = 0
        raw1, raw2, n = next_block()
        if n == 0:
            raise InputError("the input holds no reads")

        # Offsets, read lengths and an early quality gate from the first block.
        sample1 = _first_seqs(raw1, args.detect_pairs)
        sample2 = _first_seqs(raw2, args.detect_pairs)
        offsets, rates = detect_barcode_offsets(sample2, bc_map)
        r1_mode = Counter(map(len, sample1)).most_common(1)[0][0]
        r2_mode = Counter(map(len, sample2)).most_common(1)[0][0]
        sample_valid = sum(
            1 for s in sample2
            if all(s[o:o + BC_LEN] in bc_map for o in offsets)) / len(sample2)
        log(f"offsets {offsets} (match rates "
            f"{', '.join(f'{r:.3f}' for r in rates.values())}); R1 length {r1_mode}, "
            f"R2 length {r2_mode}; {sample_valid:.1%} of the first {len(sample2)} pairs "
            f"carry a valid barcode")
        if sample_valid < args.min_valid_fraction:
            raise InputError(
                f"only {sample_valid:.1%} of the first {len(sample2)} pairs carry a valid "
                f"barcode at offsets {offsets} (minimum {args.min_valid_fraction:.0%}); "
                "is this stLFR data with the barcodes still on read 2?")
        if offsets[0] < 1:
            raise InputError(f"detected barcode offsets {offsets} leave no insert in read 2")

        distinct = DistinctCounter(radix)
        totals = {"pairs": 0, "valid": 0, "kept": 0}
        pending = deque()
        max_pending = 2 * n_workers if n_workers > 1 else 1
        next_report = args.progress_pairs

        def drain_one():
            nonlocal next_report
            n_pairs, n_valid, n_kept, gz1, gz2, keys, n_kept_keys = collect(pending.popleft())
            totals["pairs"] += n_pairs
            totals["valid"] += n_valid
            totals["kept"] += n_kept
            writer.write(gz1, gz2, n_kept)
            distinct.add(keys, n_kept_keys)
            if totals["pairs"] >= next_report:
                elapsed = time.time() - t_start
                log(f"{totals['pairs']:,} pairs, {totals['kept']:,} kept, "
                    f"{totals['pairs'] / elapsed:,.0f} pairs/s")
                next_report += args.progress_pairs

        while n:
            pending.append(submit((pairs_read, offsets, raw1, raw2)))
            pairs_read += n
            del raw1, raw2
            while len(pending) >= max_pending:
                drain_one()
            raw1, raw2, n = next_block()
        while pending:
            drain_one()
        for source in sources:
            source.finish()
        writer.close()
        if pool is not None:
            pool.close()
            pool.join()
            pool = None
    finally:
        for source in sources:
            source.kill()
        if pool is not None:
            pool.terminate()
            pool.join()
        writer.close()

    elapsed = time.time() - t_start
    total = totals["pairs"]
    valid_fraction = totals["valid"] / total if total else 0.0
    stats = {
        "total_pairs": total,
        "valid_barcode_pairs": totals["valid"],
        "valid_fraction": valid_fraction,
        "kept_pairs": totals["kept"],
        "dropped_invalid": total - totals["valid"],
        "dropped_by_subsample": totals["valid"] - totals["kept"],
        "distinct_barcodes_seen": distinct.seen,
        "distinct_barcodes_kept": distinct.kept,
        "barcode_offsets": offsets,
        "barcode_offset_match_rates": [rates[o] for o in offsets],
        "insert_length_r2": offsets[0],
        "n_chunks": writer.n_chunks,
        "chunk_pairs": writer.chunk_counts,
        "chunk_files": [[os.path.basename(p) for p in pair] for pair in writer.files],
        "r1_len_mode": r1_mode,
        "r2_len_mode": r2_mode,
        "len_mode_sampled_pairs": len(sample2),
        "barcode_fraction": args.barcode_fraction,
        "seed": args.seed,
        "whitelist_size": len(entries),
        "min_valid_fraction": args.min_valid_fraction,
        "threads": threads,
        "elapsed_seconds": round(elapsed, 1),
        "pairs_per_second": round(total / elapsed) if elapsed > 0 else None,
    }
    if args.stats:
        stats_dir = os.path.dirname(args.stats)
        if stats_dir:
            os.makedirs(stats_dir, exist_ok=True)
        with open(args.stats, "w") as fh:
            json.dump(stats, fh, indent=2)
            fh.write("\n")
    log(f"done: {total:,} pairs, {totals['valid']:,} with a valid barcode "
        f"({valid_fraction:.1%}), {totals['kept']:,} kept in {writer.n_chunks} chunk(s), "
        f"{distinct.kept:,} of {distinct.seen:,} barcodes kept, {elapsed:.0f} s")
    if valid_fraction < args.min_valid_fraction:
        raise InputError(
            f"only {valid_fraction:.1%} of pairs carry a valid barcode "
            f"(minimum {args.min_valid_fraction:.0%})")
    return stats


def log(message):
    print(f"[stlfr_tag_barcodes] {message}", file=sys.stderr, flush=True)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Tag stLFR read pairs with their barcode (name#A_B_C/1) and write "
                    "FASTQ chunks, as BGI's split_barcode.pl does.")
    p.add_argument("--r1", required=True, help="read 1 FASTQ, gzipped or plain")
    p.add_argument("--r2", required=True, help="read 2 FASTQ carrying the barcodes")
    p.add_argument("--whitelist", required=True,
                   help="barcode list: sequence and integer ID per line (BGI barcode.list)")
    p.add_argument("--out-prefix", required=True,
                   help="outputs are PREFIX.chunkNNNN_1.fastq.gz and _2.fastq.gz")
    p.add_argument("--chunk-pairs", type=int, default=10_000_000,
                   help="kept pairs per output chunk [%(default)s]")
    p.add_argument("--barcode-fraction", type=float, default=1.0,
                   help="fraction of barcodes to keep, chosen by hash of (A, B, C, seed) "
                        "[%(default)s]")
    p.add_argument("--seed", type=int, default=1, help="seed for the barcode subset [%(default)s]")
    p.add_argument("--stats", help="write run statistics to this JSON file")
    p.add_argument("--threads", type=int, default=0,
                   help="CPUs to use; 0 means all available [%(default)s]")
    p.add_argument("--min-valid-fraction", type=float, default=0.5,
                   help="fail if fewer pairs than this carry a valid barcode [%(default)s]")
    p.add_argument("--compress-level", type=int, default=1,
                   help="gzip level of the output chunks [%(default)s]")
    p.add_argument("--detect-pairs", type=int, default=20_000,
                   help="leading pairs used to detect barcode offsets [%(default)s]")
    p.add_argument("--block-mb", type=float, default=8,
                   help="MB of read 1 handed to a worker at a time [%(default)s]")
    p.add_argument("--progress-pairs", type=int, default=20_000_000,
                   help="log progress every this many pairs [%(default)s]")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        run(args)
    except InputError as exc:
        print(f"[stlfr_tag_barcodes] ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

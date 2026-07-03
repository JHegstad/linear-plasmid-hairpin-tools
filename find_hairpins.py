#!/usr/bin/env python3
"""
find_hairpins.py

Detect which assembler produces HAIRPIN artefacts (a contig end that folds back
as its own reverse complement) and, for FASTA, distinguish a TRUE terminal
hairpin (the inverted repeat reaches the contig end) from an INTERNAL inverted
repeat (rRNA operons etc., which are not hairpins).

Evidence differs by file type:

  GFA  (authoritative for graph-encoded folds):
       a hairpin is an opposite-strand self-link  L <seg> + <seg> -
       Per GFA we also print total links and same-strand self-links, so you can
       confirm the parser sees links (0 hairpin links is then a real result, not
       a parsing miss). NB: most long-read assemblers bake the fold into the
       contig SEQUENCE and emit it linearly, so they legitimately have 0 hairpin
       links even when the FASTA shows a terminal fold.

  FASTA (heuristic, with localisation):
       within a terminal window we find k-mers whose reverse complement also
       occurs in the window (an inverted repeat). We locate the symmetry axis
       (fold apex) and measure whether the outer arm reaches the contig terminus
       (within --edge-tol). If it does -> TERMINAL hairpin; otherwise INTERNAL.

Files are grouped by assembler (prefix before the first "_<number>").

Usage
-----
    python find_hairpins.py [--dir .] [-k 31] [--window 50000]
                            [--min-shared 25] [--edge-tol 100] [--gfa-only]

Outputs (written to --dir):
    hairpin_report.txt   human-readable, grouped by assembler
    hairpin_report.tsv   one row per detected inverted repeat
"""

import argparse
import glob
import os
import re
from collections import defaultdict, Counter

_COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def rev_comp(s: bytes) -> bytes:
    return s.translate(_COMP)[::-1]


# ── FASTA reading ─────────────────────────────────────────────────────────────

def read_fasta(path):
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name = line[1:].split()[0].decode("ascii", "replace")
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


# ── FASTA inverted-repeat localisation ────────────────────────────────────────

def analyse_end(seq: bytes, k: int, window: int, end: str, edge_tol: int,
                max_kmer_hits: int = 200):
    """
    Look for a reverse-complement inverted repeat in the terminal window at the
    given end ('5' or '3'). Locate the dominant symmetry axis and measure the
    outer arm.

    Returns None or dict:
      support  : number of supporting RC k-mer pairs at the dominant axis
      tir_span : bp spanned by the inverted repeat (both arms + loop)
      touches  : True if the outer arm reaches the contig terminus (<= edge_tol)
      gap      : bp from the terminus to the outer arm (0 = flush with the end)
    """
    L = len(seq)
    if L < k:
        return None
    if end == "3":
        off = max(0, L - window)
        w = seq[off:]
    else:  # "5"
        off = 0
        w = seq[:min(window, L)]
    n = len(w)
    if n < k:
        return None

    pos = defaultdict(list)
    for i in range(n - k + 1):
        pos[w[i:i + k]].append(i)

    # centre (axis) -> list of (i, j) supporting pairs
    centres = defaultdict(list)
    for i in range(n - k + 1):
        js = pos.get(rev_comp(w[i:i + k]))
        if not js or len(js) > max_kmer_hits:
            continue
        for j in js:
            if j < i:
                continue
            centres[(i + j + k - 1) // 2].append((i, j))
    if not centres:
        return None

    c0 = max(centres, key=lambda c: len(centres[c]))
    pairs = []
    for c in (c0 - 1, c0, c0 + 1):          # allow ±1 for rounding
        pairs += centres.get(c, [])
    positions = set()
    for i, j in pairs:
        positions.add(i)
        positions.add(j)
    support = len(pairs)
    minp, maxp = min(positions), max(positions)
    tir_span = (maxp - minp) + k

    lefts = [i for i, j in pairs]
    rights = [j for i, j in pairs]
    # arm1 = 5'-side arm, arm2 = 3'-side arm (global 0-based half-open)
    arm1 = (off + min(lefts), off + max(lefts) + k)
    arm2 = (off + min(rights), off + max(rights) + k)

    outer_global = off + maxp + k          # 3' edge of outermost k-mer
    inner_global = off + minp              # 5' edge of innermost k-mer
    if end == "3":
        gap = max(0, L - outer_global)
        touches = gap <= edge_tol
    else:
        gap = max(0, inner_global)
        touches = gap <= edge_tol
    return {"support": support, "tir_span": tir_span,
            "touches": touches, "gap": gap,
            "arm1": arm1, "arm2": arm2, "axis": off + c0}


def scan_fasta(path, k, window, min_shared, edge_tol, keep_seq=False):
    n_contigs = 0
    hits = []          # dicts with per-contig detail
    for name, seq in read_fasta(path):
        n_contigs += 1
        seq = seq.upper()
        best = None
        for end in ("5", "3"):
            r = analyse_end(seq, k, window, end, edge_tol)
            if r is None or r["support"] < min_shared:
                continue
            # prefer an end that TOUCHES the terminus, then higher support
            key = (r["touches"], r["support"])
            if best is None or key > best[0]:
                best = (key, end, r)
        if best is not None:
            _, end, r = best
            hit = {"contig": name, "length": len(seq), "end": end,
                   "tir_span": r["tir_span"], "gap": r["gap"],
                   "support": r["support"], "terminal": r["touches"],
                   "arm1": r["arm1"], "arm2": r["arm2"], "axis": r["axis"]}
            if keep_seq:
                hit["seq"] = seq
            hits.append(hit)
    return n_contigs, hits


# ── GFA ───────────────────────────────────────────────────────────────────────

def scan_gfa(path):
    """Return dict: n_seg, n_links, n_self_same, n_self_opp, hairpin_segs."""
    segs = set()
    n_links = 0
    self_same = 0
    hairpin_segs = []
    with open(path) as fh:
        for line in fh:
            if not line:
                continue
            t = line[0]
            if t == "S":
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    segs.add(f[1])
            elif t == "L":
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                n_links += 1
                if f[1] == f[3]:
                    if f[2] == f[4]:
                        self_same += 1
                    else:
                        hairpin_segs.append(f[1])
    return {"n_seg": len(segs), "n_links": n_links, "self_same": self_same,
            "self_opp": len(hairpin_segs),
            "hairpin_segs": sorted(set(hairpin_segs))}


# ── grouping ──────────────────────────────────────────────────────────────────

def assembler_of(filename: str) -> str:
    base = os.path.basename(filename)
    m = re.match(r"([A-Za-z][A-Za-z0-9]*?)_\d+\.", base)
    return m.group(1) if m else base.split("_")[0].split(".")[0]


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Detect hairpins per assembler; GFA links authoritative, "
                    "FASTA terminal self-RC localised (terminal vs internal).")
    ap.add_argument("--dir", default=".", help="Directory to scan (default: .)")
    ap.add_argument("-k", type=int, default=31, help="k-mer size (FASTA)")
    ap.add_argument("--window", type=int, default=50000,
                    help="Terminal window (bp) scanned at each contig end")
    ap.add_argument("--min-shared", type=int, default=25,
                    help="Min supporting RC k-mer pairs to report an IR")
    ap.add_argument("--edge-tol", type=int, default=100,
                    help="Max gap (bp) from terminus to call a fold TERMINAL")
    ap.add_argument("--gfa-only", action="store_true",
                    help="Only use GFA files (skip the FASTA heuristic)")
    ap.add_argument("--extract-hairpins", action="store_true",
                    help="Write flagged TERMINAL-hairpin contigs to "
                         "hairpin_contigs.fasta (+ hairpin_tir_coords.tsv) "
                         "for dotplot/alignment.")
    args = ap.parse_args()

    gfas = sorted(glob.glob(os.path.join(args.dir, "*.gfa")))
    fastas = sorted(glob.glob(os.path.join(args.dir, "*.fasta")))
    gfa_stems = {os.path.splitext(os.path.basename(g))[0] for g in gfas}

    report = []

    def emit(s=""):
        print(s)
        report.append(s)

    emit(f"Scanning {args.dir}")
    emit(f"  {len(gfas)} GFA file(s); {len(fastas)} FASTA file(s)")
    emit(f"  FASTA params: k={args.k}, window={args.window}, "
         f"min_shared={args.min_shared}, edge_tol={args.edge_tol}")
    emit()

    file_rows = []     # per-file summary rows
    ir_rows = []       # per-inverted-repeat rows (for TSV)

    # ── GFA diagnostics ───────────────────────────────────────────────────────
    emit("GFA link diagnostics (confirms the parser sees links)")
    emit("-" * 72)
    emit(f"{'file':<22}{'segs':>6}{'links':>7}{'self=':>7}{'hairpin':>9}  segs")
    for g in gfas:
        s = scan_gfa(g)
        emit(f"{os.path.basename(g):<22}{s['n_seg']:>6}{s['n_links']:>7}"
             f"{s['self_same']:>7}{s['self_opp']:>9}  "
             f"{','.join(s['hairpin_segs']) if s['hairpin_segs'] else '-'}")
        file_rows.append(dict(file=os.path.basename(g), assembler=assembler_of(g),
                              type="gfa", n_units=s["n_seg"],
                              n_hairpin=s["self_opp"], n_internal=0,
                              gfa=True))
        for seg in s["hairpin_segs"]:
            ir_rows.append([os.path.basename(g), assembler_of(g), "gfa", seg,
                            "", "hairpin-link", "", "", ""])
    emit()

    # ── FASTA localised scan ──────────────────────────────────────────────────
    if not args.gfa_only:
        emit("FASTA terminal-fold scan (TERMINAL = reaches contig end)")
        emit("-" * 72)
        emit(f"{'file':<22}{'ctgs':>5}{'term':>5}{'int':>5}  "
             f"terminal-hairpin contigs (len, end, tir_bp, gap_bp, support)")
        extracted = []          # (stem, hit) for terminal hits, if extracting
        for fa in fastas:
            n_contigs, hits = scan_fasta(fa, args.k, args.window,
                                         args.min_shared, args.edge_tol,
                                         keep_seq=args.extract_hairpins)
            term = [h for h in hits if h["terminal"]]
            if args.extract_hairpins:
                stem = os.path.splitext(os.path.basename(fa))[0]
                for h in term:
                    extracted.append((stem, h))
            intern = [h for h in hits if not h["terminal"]]
            has_gfa = os.path.splitext(os.path.basename(fa))[0] in gfa_stems
            det = "; ".join(
                f"{h['contig']}({h['length']}bp,{h['end']}',"
                f"{h['tir_span']}bp,gap{h['gap']},{h['support']}rc)"
                for h in term) or "-"
            flag = "  <== HAIRPIN" if term else ""
            emit(f"{os.path.basename(fa):<22}{n_contigs:>5}{len(term):>5}"
                 f"{len(intern):>5}  {det}{flag}")
            file_rows.append(dict(file=os.path.basename(fa),
                                  assembler=assembler_of(fa),
                                  type="fasta*" if has_gfa else "fasta",
                                  n_units=n_contigs, n_hairpin=len(term),
                                  n_internal=len(intern), gfa=False))
            for h in hits:
                ir_rows.append([os.path.basename(fa), assembler_of(fa),
                                "fasta", h["contig"], h["length"],
                                "terminal" if h["terminal"] else "internal",
                                f"{h['end']}'", h["tir_span"], h["gap"]])
        emit()

    # ── per-assembler verdict ─────────────────────────────────────────────────
    by_asm = defaultdict(lambda: {"files": 0, "hp": 0, "gfa": False})
    for r in file_rows:
        a = by_asm[r["assembler"]]
        a["files"] += 1
        if r["n_hairpin"]:
            a["hp"] += 1
        if r["gfa"]:
            a["gfa"] = True

    emit("Per-assembler verdict (terminal hairpin present?)")
    emit("-" * 72)
    emit(f"{'assembler':<14}{'files':>6}{'w/hairpin':>11}  best evidence")
    for asm in sorted(by_asm):
        a = by_asm[asm]
        ev = "GFA self-link (authoritative)" if a["gfa"] else "FASTA terminal fold"
        verdict = "HAIRPIN" if a["hp"] else "clean"
        emit(f"{asm:<14}{a['files']:>6}{a['hp']:>11}  {verdict:<8} [{ev}]")
    emit()
    emit("Notes:")
    emit("  * GFA 'hairpin' counts opposite-strand self-links (L n + n -).")
    emit("    0 there with folds in the FASTA just means the assembler wrote the")
    emit("    fold as sequence, not as a graph edge.")
    emit("  * FASTA 'term' = fold reaches the contig end (gap <= edge_tol) = a")
    emit("    real terminal hairpin; 'int' = internal inverted repeat (e.g. rRNA).")
    emit("  * A feature reproduced on the same contig by many assemblers is")
    emit("    likely genuine biology (linear-plasmid terminal inverted repeats).")

    # ── extract flagged terminal-hairpin contigs ──────────────────────────────
    if args.extract_hairpins and not args.gfa_only:
        fa_out = os.path.join(args.dir, "hairpin_contigs.fasta")
        co_out = os.path.join(args.dir, "hairpin_tir_coords.tsv")
        n_written = 0
        with open(fa_out, "w") as fo, open(co_out, "w") as co:
            co.write("file\tcontig\tlength\tend\ttir_span\tgap\taxis\t"
                     "arm1_start\tarm1_end\tarm2_start\tarm2_end\tsupport\n")
            for stem, h in extracted:
                a1s, a1e = h["arm1"]
                a2s, a2e = h["arm2"]
                name = f"{stem}__{h['contig']}"
                header = (f">{name} len={h['length']} end={h['end']}prime "
                          f"tir_span={h['tir_span']} gap={h['gap']} "
                          f"arm1={a1s + 1}-{a1e} arm2={a2s + 1}-{a2e} "
                          f"support={h['support']}")
                seq = h["seq"]
                fo.write(header + "\n")
                for i in range(0, len(seq), 80):
                    fo.write(seq[i:i + 80].decode("ascii", "replace") + "\n")
                # 1-based inclusive coords in the TSV
                co.write(f"{stem}\t{h['contig']}\t{h['length']}\t{h['end']}\t"
                         f"{h['tir_span']}\t{h['gap']}\t{h['axis']}\t"
                         f"{a1s + 1}\t{a1e}\t{a2s + 1}\t{a2e}\t{h['support']}\n")
                n_written += 1
        emit(f"Extracted {n_written} terminal-hairpin contig(s) -> "
             f"{os.path.basename(fa_out)} (+ {os.path.basename(co_out)})")
        emit()

    # ── write outputs ─────────────────────────────────────────────────────────
    txt = os.path.join(args.dir, "hairpin_report.txt")
    with open(txt, "w") as fh:
        fh.write("\n".join(report) + "\n")
    tsv = os.path.join(args.dir, "hairpin_report.tsv")
    with open(tsv, "w") as fh:
        fh.write("file\tassembler\tsource\tunit\tlength\tclass\tend\t"
                 "tir_span\tgap\n")
        for row in ir_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nReport written to : {txt}")
    print(f"Table written to  : {tsv}")


if __name__ == "__main__":
    main()

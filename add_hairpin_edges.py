#!/usr/bin/env python3
"""
add_hairpin_edges.py

Detect terminal fold-backs (hairpins) in a raw assembler GFA (e.g. plassembler)
and add the SAME hairpin link Autocycler uses to represent them:

    3' fold-back  ->  L <seg> + <seg> - 0M     (forward end folds to own reverse)
    5' fold-back  ->  L <seg> - <seg> + 0M     (reverse end folds to own forward)

This matches Autocycler's convention (unitig.rs::hairpin_end / hairpin_start,
written by unitig_graph.rs::get_links_for_gfa), so the hairpin becomes visible
to Bandage, to Autocycler's topology(), and to find_hairpins/the classifier.

SCOPE / honesty
---------------
This ANNOTATES the graph — it adds the link but does NOT alter segment
sequences. Autocycler's internal representation additionally collapses the
duplicated fold-back arm into a single unitig (so the segment holds each arm
once and the 0M link expresses the turn). Faithfully reproducing that from a
raw contig requires splitting the contig at the fold apex (graph surgery); that
is deliberately NOT done here. The added link is a correct topological
annotation; the sequence still contains both arms of the fold.

A hairpin is detected exactly like find_hairpins: within a terminal window we
look for k-mers whose reverse complement is also present (an inverted repeat),
locate the symmetry axis, and require the outer arm to reach the contig end
(gap <= --edge-tol).

Usage
-----
    python add_hairpin_edges.py plassembler_01.gfa [more.gfa ...]
        [-k 31] [--window 50000] [--min-shared 25] [--edge-tol 100]
        [--overlap] [--suffix .hairpins] [--inplace]

Writes <stem><suffix>.gfa next to each input (default suffix '.hairpins'),
unless --inplace is given.
"""

import argparse
import os
import sys
from collections import defaultdict

_COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def rev_comp(s: bytes) -> bytes:
    return s.translate(_COMP)[::-1]


def analyse_end(seq: bytes, k: int, window: int, end: str, edge_tol: int,
                max_kmer_hits: int = 200):
    """Return dict(support, arm_len, touches, gap) for a terminal fold, or None."""
    L = len(seq)
    if L < k:
        return None
    off = max(0, L - window) if end == "3" else 0
    w = seq[off:] if end == "3" else seq[:min(window, L)]
    n = len(w)
    if n < k:
        return None
    pos = defaultdict(list)
    for i in range(n - k + 1):
        pos[w[i:i + k]].append(i)
    centres = defaultdict(list)
    for i in range(n - k + 1):
        js = pos.get(rev_comp(w[i:i + k]))
        if not js or len(js) > max_kmer_hits:
            continue
        for j in js:
            if j >= i:
                centres[(i + j + k - 1) // 2].append((i, j))
    if not centres:
        return None
    c0 = max(centres, key=lambda c: len(centres[c]))
    pairs = centres.get(c0 - 1, []) + centres.get(c0, []) + centres.get(c0 + 1, [])
    positions = {p for pr in pairs for p in pr}
    minp, maxp = min(positions), max(positions)
    lefts = [i for i, _ in pairs]
    rights = [j for _, j in pairs]
    arm_len = min((max(lefts) - min(lefts)) + k, (max(rights) - min(rights)) + k)
    outer = off + maxp + k
    inner = off + minp
    if end == "3":
        gap = max(0, L - outer)
    else:
        gap = max(0, inner)
    return {"support": len(pairs), "arm_len": arm_len,
            "touches": gap <= edge_tol, "gap": gap}


def assembler_of(name):
    import re
    m = re.match(r"([A-Za-z][A-Za-z0-9]*?)_\d+\.", os.path.basename(name))
    return m.group(1) if m else os.path.basename(name).split("_")[0].split(".")[0]


def process_gfa(path, args, records, summaries):
    lines = []
    seqs = {}
    existing_links = set()
    with open(path) as fh:
        for line in fh:
            lines.append(line.rstrip("\n"))
            if line[:1] == "S":
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[2] != "*":
                    seqs[f[1]] = f[2].encode("ascii", "replace").upper()
            elif line[:1] == "L":
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5:
                    existing_links.add((f[1], f[2], f[3], f[4]))

    n_links_in = len(existing_links)
    new_links = []          # (seg, from_strand, to_strand, overlap, info)
    for seg, seq in seqs.items():
        for end in ("5", "3"):
            r = analyse_end(seq, args.k, args.window, end, args.edge_tol)
            if r is None or r["support"] < args.min_shared or not r["touches"]:
                continue
            # 3' fold -> L seg + seg - ; 5' fold -> L seg - seg +
            sa, sb = ("+", "-") if end == "3" else ("-", "+")
            if (seg, sa, seg, sb) in existing_links:
                continue
            ov = f"{r['arm_len']}M" if args.overlap else "0M"
            new_links.append((seg, sa, sb, ov,
                              f"end={end}' arm={r['arm_len']}bp gap={r['gap']} "
                              f"support={r['support']}"))
            records.append({"file": os.path.basename(path), "segment": seg,
                            "from_strand": sa, "to_strand": sb, "overlap": ov,
                            "end": f"{end}'", "arm_len": r["arm_len"],
                            "gap": r["gap"], "support": r["support"],
                            "seg_len": len(seq)})
            existing_links.add((seg, sa, seg, sb))

    hp_segs = sorted({seg for seg, _sa, _sb, _ov, _i in new_links})
    max_arm = max((r["arm_len"] for r in records
                   if r["file"] == os.path.basename(path)), default=0)

    if not new_links:
        out = ""
        print(f"{os.path.basename(path)}: no terminal hairpins detected — "
              f"no links added")
    else:
        out = path if args.inplace else (
            os.path.splitext(path)[0] + args.suffix + ".gfa")
        with open(out, "w") as fh:
            for ln in lines:
                fh.write(ln + "\n")
            for seg, sa, sb, ov, _info in new_links:
                fh.write(f"L\t{seg}\t{sa}\t{seg}\t{sb}\t{ov}\n")
        print(f"{os.path.basename(path)}: added {len(new_links)} hairpin "
              f"link(s) -> {os.path.basename(out)}")
        for seg, sa, sb, ov, info in new_links:
            print(f"    L {seg} {sa} {seg} {sb} {ov}   [{info}]")

    summaries.append({"file": os.path.basename(path),
                      "assembler": assembler_of(path),
                      "n_segments": len(seqs), "n_links_in": n_links_in,
                      "n_hairpins_added": len(new_links),
                      "max_arm_len": max_arm,
                      "hairpin_segments": ",".join(hp_segs) if hp_segs else "-",
                      "out_gfa": os.path.basename(out) if out else "-"})
    return out, len(new_links)


def main():
    ap = argparse.ArgumentParser(
        description="Add Autocycler-style hairpin links to raw assembler GFAs.")
    ap.add_argument("gfa", nargs="+", help="Input GFA file(s)")
    ap.add_argument("-k", type=int, default=31, help="k-mer size")
    ap.add_argument("--window", type=int, default=50000,
                    help="Terminal window (bp) scanned at each segment end")
    ap.add_argument("--min-shared", type=int, default=25,
                    help="Min supporting RC k-mer pairs to call a hairpin")
    ap.add_argument("--edge-tol", type=int, default=100,
                    help="Max gap (bp) from terminus to accept the fold")
    ap.add_argument("--overlap", action="store_true",
                    help="Emit the measured arm length as the link overlap "
                         "(e.g. 3000M) instead of Autocycler's 0M")
    ap.add_argument("--suffix", default=".hairpins",
                    help="Suffix for output files (default: .hairpins)")
    ap.add_argument("--inplace", action="store_true",
                    help="Overwrite the input GFA instead of writing a copy")
    ap.add_argument("--report", default=None,
                    help="Path for the per-link report TSV "
                         "(default: hairpin_edges.tsv next to the first GFA)")
    ap.add_argument("--summary", default=None,
                    help="Path for the per-GFA summary TSV "
                         "(default: hairpin_summary.tsv next to the first GFA)")
    ap.add_argument("--append", action="store_true",
                    help="Append to existing report/summary TSVs instead of "
                         "overwriting them (header written only for a new file)")
    args = ap.parse_args()

    total = 0
    records = []          # one row per added link
    summaries = []        # one row per GFA
    for path in args.gfa:
        if not os.path.exists(path):
            print(f"SKIP (not found): {path}", file=sys.stderr)
            continue
        _, n = process_gfa(path, args, records, summaries)
        total += n

    base = os.path.dirname(os.path.abspath(args.gfa[0]))
    report = args.report or os.path.join(base, "hairpin_edges.tsv")
    summary = args.summary or os.path.join(base, "hairpin_summary.tsv")
    mode = "a" if args.append else "w"

    def _need_header(p):
        # write header unless appending to a non-empty existing file
        return not (args.append and os.path.exists(p) and os.path.getsize(p) > 0)

    # ── per-link detail report ────────────────────────────────────────────────
    dcols = ["file", "segment", "seg_len", "from_strand", "to_strand",
             "overlap", "end", "arm_len", "gap", "support"]
    with open(report, mode) as fh:
        if _need_header(report):
            fh.write("# Autocycler-style hairpin links (add_hairpin_edges.py)\n")
            fh.write(f"# k={args.k} window={args.window} "
                     f"min_shared={args.min_shared} edge_tol={args.edge_tol} "
                     f"overlap={'arm_len' if args.overlap else '0M'}\n")
            fh.write("\t".join(dcols) + "\n")
        for r in records:
            fh.write("\t".join(str(r[c]) for c in dcols) + "\n")

    # ── per-GFA summary (all GFAs, incl. those with 0 hairpins) ────────────────
    scols = ["file", "assembler", "n_segments", "n_links_in",
             "n_hairpins_added", "max_arm_len", "hairpin_segments", "out_gfa"]
    with open(summary, mode) as fh:
        if _need_header(summary):
            fh.write("\t".join(scols) + "\n")
        for s in summaries:
            fh.write("\t".join(str(s[c]) for c in scols) + "\n")

    verb = "appended to" if args.append else "written to"
    print(f"\nDone. {total} hairpin link(s) added across {len(summaries)} "
          f"GFA file(s).")
    print(f"Per-link report {verb} : {report}")
    print(f"Per-GFA summary {verb} : {summary}")


if __name__ == "__main__":
    main()

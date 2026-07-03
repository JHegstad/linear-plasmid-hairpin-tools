#!/usr/bin/env python3
"""
autocycler_dotplot_classify.py

Classify Autocycler clusters / replicons as CIRCULAR or LINEAR by reading the
GFA graph directly (NOT by analysing dotplot PNG pixels).

Two input modes
---------------
1. Per-cluster (default):
     <search>/*_autocycler/autocycler_out/clustering/qc_pass/cluster_*/<gfa>
   where <gfa> is chosen with --gfa-name (default 1_untrimmed.gfa; you can
   point it at 5_final.gfa for the trimmed/resolved graph).
   Each `P` line is one input assembly's contig, given as an ordered unitig
   path. A contig is CIRCULAR iff its path closes into a loop (a link from the
   path end back to its start). The cluster topology is the majority vote.

2. Consensus (--consensus):
     <search>/*_autocycler/autocycler_out/consensus_assembly.gfa
   The final combined assembly. It has no `P` lines, so each replicon is a
   connected component of unitigs. Topology is decided from the graph
   structure (Autocycler's own rule): a component that forms a simple
   circular loop is CIRCULAR; otherwise LINEAR (open or hairpin), or
   FRAGMENTED if it is several unitigs that do not form a clean loop.

Why not classify from the PNG?
------------------------------
The old version counted red vs blue dotplot pixels. Red = reverse-complement
k-mer matches (dotplot.rs REVERSE_DOT_COLOUR), which reflects strand
disagreement between assemblies and inverted repeats/hairpins — NOT topology.
A circular plasmid whose assemblies are in opposite orientations shows lots of
red (looked "linear"); a linear plasmid with a short/absent hairpin shows
little red (looked "circular"). Plus the blue self-vs-self diagonals dominate
the counts, so any fixed threshold is arbitrary.

GFA topology encoding (unitig_graph.rs)
---------------------------------------
  H  VN:Z:1.0  KM:i:<k>
  S  <id>  <seq>  <tags...>
  L  <a> <a_strand> <b> <b_strand>  0M
  P  <seq_id> <u1><s1>,<u2><s2>,...  *  LN:i:.. FN:Z:<file> HD:Z:<header>

  * same-strand self/closing link (L n + n + / L n - n -)  -> circularising
  * opposite-strand self link     (L n + n - / L n - n +)  -> hairpin (linear)

`--strict` (per-cluster mode) additionally computes the graph-structure
topology of the whole cluster graph (mirroring Autocycler's
`component_is_circular_loop`) and flags any cluster where the P-line vote and
the graph structure disagree, so you can eyeball it.

Usage
-----
    python autocycler_dotplot_classify.py [--search-dir DIR]
                                          [--gfa-name 1_untrimmed.gfa]
                                          [--strict]
                                          [--consensus]
                                          [--no-dotplot | --make-dotplot-only]
"""

import argparse
import glob
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


# ── GFA parsing ───────────────────────────────────────────────────────────────

def parse_gfa(gfa_path: Path):
    """
    Parse an Autocycler GFA.

    Returns (segments, links, paths):
      segments : dict seg_id -> length (bp; from sequence, else LN:i tag, else 0)
      links    : set of (from_seg, from_strand, to_seg, to_strand)
                 Autocycler writes each link and its RC mirror, so this set is
                 orientation-complete.
      paths    : list of {'id','path':[(seg,strand),...],'filename','header','length'}
    """
    segments = {}
    links = set()
    paths = []
    with open(gfa_path) as fh:
        for line in fh:
            if not line:
                continue
            t = line[0]
            if t not in "SLP":
                continue
            f = line.rstrip("\n").split("\t")
            if t == "S" and len(f) >= 3:
                seg = f[1]
                seq = f[2]
                length = 0
                if seq and seq != "*":
                    length = len(seq)
                else:
                    for tag in f[3:]:
                        if tag.startswith("LN:i:"):
                            try:
                                length = int(tag[5:])
                            except ValueError:
                                pass
                segments[seg] = length
            elif t == "L" and len(f) >= 5:
                links.add((f[1], f[2], f[3], f[4]))
            elif t == "P" and len(f) >= 3:
                seg_field = f[2].strip()
                if not seg_field or seg_field == "*":
                    continue
                path, ok = [], True
                for tok in seg_field.split(","):
                    tok = tok.strip()
                    if len(tok) < 2 or tok[-1] not in "+-":
                        ok = False
                        break
                    path.append((tok[:-1], tok[-1]))
                if not ok or not path:
                    continue
                rec = {"id": f[1], "path": path, "filename": "",
                       "header": "", "length": None}
                for tag in f[3:]:
                    if tag.startswith("FN:Z:"):
                        rec["filename"] = tag[5:]
                    elif tag.startswith("HD:Z:"):
                        rec["header"] = tag[5:]
                    elif tag.startswith("LN:i:"):
                        try:
                            rec["length"] = int(tag[5:])
                        except ValueError:
                            pass
                paths.append(rec)
    # Make sure every segment referenced only in links/paths still exists.
    for (a, _, b, _) in links:
        segments.setdefault(a, 0)
        segments.setdefault(b, 0)
    return segments, links, paths


# ── P-line (input-assembly) circularity ───────────────────────────────────────

def contig_is_circular(path, links) -> bool:
    """Circular iff a link closes the path end back to its start (same travel)."""
    if not path:
        return False
    first_seg, first_strand = path[0]
    last_seg, last_strand = path[-1]
    return (last_seg, last_strand, first_seg, first_strand) in links


def has_hairpin_link(links) -> bool:
    return any(a == b and sa != sb for (a, sa, b, sb) in links)


# ── graph-structure topology (mirrors UnitigGraph::component_is_circular_loop) ─

def _build_adjacency(links):
    """Return per-unitig link tables needed for the circular-loop walk."""
    fnext = defaultdict(list)   # leaving u on '+'  -> list of (to, to_strand)
    rnext = defaultdict(list)   # leaving u on '-'  -> list of (to, to_strand)
    fprev = defaultdict(int)    # links arriving at u on '+'
    rprev = defaultdict(int)    # links arriving at u on '-'
    undirected = defaultdict(set)
    for (a, sa, b, sb) in links:
        (fnext if sa == "+" else rnext)[a].append((b, sb))
        if sb == "+":
            fprev[b] += 1
        else:
            rprev[b] += 1
        undirected[a].add(b)
        undirected[b].add(a)
    return fnext, rnext, fprev, rprev, undirected


def connected_components(segments, undirected):
    visited, comps = set(), []
    for seg in segments:
        if seg in visited:
            continue
        stack, comp = [seg], []
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            for nb in undirected.get(cur, ()):
                if nb not in visited:
                    stack.append(nb)
        comps.append(sorted(comp))
    return sorted(comps)


def component_is_circular_loop(component, fnext, rnext, fprev, rprev) -> bool:
    """
    Faithful port of Autocycler's component_is_circular_loop: start at the
    lowest-numbered unitig on the forward strand and walk forward links; every
    unitig on the loop must have exactly one link on each of its four sides,
    and the walk must return to the start after covering the whole component.
    """
    if not component:
        return False
    first = component[0]
    num, strand = first, "+"
    visited = set()
    while num != first or not visited:
        if num in visited:
            return False
        visited.add(num)
        if (len(fnext.get(num, [])) != 1 or len(rnext.get(num, [])) != 1 or
                fprev.get(num, 0) != 1 or rprev.get(num, 0) != 1):
            return False
        nxt = fnext[num][0] if strand == "+" else rnext[num][0]
        num, strand = nxt
    return len(visited) == len(component)


def classify_components(segments, links):
    """
    Classify every connected component of the graph.
    Returns list of dicts: {'segs','length','topology','hairpin'}
    topology in {'circular','linear','fragmented'}.
    """
    fnext, rnext, fprev, rprev, undirected = _build_adjacency(links)
    out = []
    for comp in connected_components(segments, undirected):
        cset = set(comp)
        length = sum(segments.get(s, 0) for s in comp)
        hp = any(a == b and sa != sb for (a, sa, b, sb) in links if a in cset)
        touches = any((a in cset or b in cset) for (a, _, b, _) in links)
        if component_is_circular_loop(comp, fnext, rnext, fprev, rprev):
            topo = "circular"
        elif len(comp) == 1 and not touches:
            topo = "linear"          # isolated unitig, open ends
        elif len(comp) == 1:
            topo = "linear"          # single unitig with a hairpin link
        else:
            topo = "fragmented"      # several unitigs, not a clean loop
        out.append({"segs": comp, "length": length, "topology": topo,
                    "hairpin": hp})
    return out


def graph_topology_label(components):
    """Collapse component classifications into one label for a cluster graph."""
    if not components:
        return "unknown"
    # Use the largest component (the replicon) as representative.
    main = max(components, key=lambda c: (c["length"], len(c["segs"])))
    return main["topology"]


# ── per-cluster (P-line) classification ───────────────────────────────────────

def classify_cluster(gfa_path: Path, strict: bool = False):
    try:
        segments, links, paths = parse_gfa(gfa_path)
    except OSError as e:
        return {"topology": "unknown", "n_seqs": 0, "n_circular": 0,
                "n_linear": 0, "hairpin": False, "graph_topology": "",
                "review": False, "detail": f"read error: {e}"}

    graph_topo = graph_topology_label(classify_components(segments, links))

    if not paths:
        # No input-assembly paths (e.g. a resolved single-contig graph).
        # Fall back to the graph structure.
        topo = graph_topo if graph_topo != "unknown" else "unknown"
        return {"topology": topo, "n_seqs": 0, "n_circular": 0, "n_linear": 0,
                "hairpin": has_hairpin_link(links), "graph_topology": graph_topo,
                "review": False,
                "detail": f"no P lines; graph structure = {graph_topo}"}

    n_seqs = len(paths)
    n_circular = sum(1 for p in paths if contig_is_circular(p["path"], links))
    n_linear = n_seqs - n_circular
    hairpin = has_hairpin_link(links)

    if n_circular == n_seqs:
        topology = "circular"
    elif n_circular == 0:
        topology = "linear"
    else:
        topology = "mixed"

    detail = f"{n_circular}/{n_seqs} contigs circular"
    if topology == "linear" and hairpin:
        detail += "; hairpin link present"

    review = False
    if strict:
        vote_circular = topology == "circular"
        struct_circular = graph_topo == "circular"
        if graph_topo in ("circular", "linear", "fragmented") and \
                vote_circular != struct_circular:
            review = True
            detail += f"; STRICT: graph structure = {graph_topo} (disagrees)"
        else:
            detail += f"; graph structure = {graph_topo}"

    return {"topology": topology, "n_seqs": n_seqs, "n_circular": n_circular,
            "n_linear": n_linear, "hairpin": hairpin,
            "graph_topology": graph_topo, "review": review, "detail": detail}


# ── autocycler dotplot runner (optional, visual QC only) ──────────────────────

def run_dotplot(gfa: Path, png: Path) -> bool:
    cmd = ["autocycler", "dotplot", "-i", str(gfa), "-o", str(png)]
    print(f"    $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    ERROR (exit {result.returncode}):", file=sys.stderr)
        if result.stderr:
            print(f"    {result.stderr.strip()}", file=sys.stderr)
        return False
    return True


# ── discovery ─────────────────────────────────────────────────────────────────

def find_qc_dirs(search_dir: Path):
    pat = str(search_dir / "**" / "*_autocycler" / "autocycler_out" /
              "clustering" / "qc_pass")
    dirs = sorted(Path(p) for p in glob.glob(pat, recursive=True))
    if not dirs:
        pat = str(search_dir / "*_autocycler" / "autocycler_out" /
                  "clustering" / "qc_pass")
        dirs = sorted(Path(p) for p in glob.glob(pat))
    return dirs


def find_consensus_gfas(search_dir: Path):
    pat = str(search_dir / "**" / "*_autocycler" / "autocycler_out" /
              "consensus_assembly.gfa")
    gfas = sorted(Path(p) for p in glob.glob(pat, recursive=True))
    if not gfas:
        pat = str(search_dir / "*_autocycler" / "autocycler_out" /
                  "consensus_assembly.gfa")
        gfas = sorted(Path(p) for p in glob.glob(pat))
    # Also accept the search dir pointing straight at an autocycler_out folder.
    direct = search_dir / "consensus_assembly.gfa"
    if direct.exists() and direct not in gfas:
        gfas.append(direct)
    return gfas


def sample_name(path: Path) -> str:
    for part in path.parts:
        if part.endswith("_autocycler"):
            return part[: -len("_autocycler")]
    return path.parts[-2] if len(path.parts) >= 2 else str(path)


# ── consensus mode ────────────────────────────────────────────────────────────

def run_consensus(args) -> None:
    gfas = find_consensus_gfas(args.search_dir)
    if not gfas:
        print(f"No 'consensus_assembly.gfa' found under '{args.search_dir}'.",
              file=sys.stderr)
        sys.exit(1)

    report = []                       # collect the human-readable output

    def emit(line=""):                # print to console AND save to report
        print(line)
        report.append(line)

    emit(f"Found {len(gfas)} consensus_assembly.gfa file(s).\n")
    results = []
    for gfa in gfas:
        sample = sample_name(gfa)
        emit("─" * 60)
        emit(f"{gfa}  (sample: {sample})")
        segments, links, _ = parse_gfa(gfa)
        comps = classify_components(segments, links)
        comps.sort(key=lambda c: -c["length"])
        for i, c in enumerate(comps, 1):
            hp = " (hairpin)" if c["hairpin"] and c["topology"] != "circular" else ""
            emit(f"  replicon {i:2d}: {c['topology'].upper():10s} "
                 f"{c['length']:>10,} bp  unitigs={len(c['segs'])}{hp}")
            results.append(dict(sample=sample, replicon=i,
                                topology=c["topology"], length=c["length"],
                                unitigs=len(c["segs"]), hairpin=c["hairpin"],
                                gfa=str(gfa)))

    emit(f"\n{'='*60}\nSUMMARY (consensus replicons)\n{'='*60}")
    counts = Counter(r["topology"] for r in results)
    emit("Totals: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    # Human-readable report (mirrors the console output).
    txt = args.search_dir / "consensus_topology.txt"
    with open(txt, "w") as fh:
        fh.write("\n".join(report) + "\n")

    # Machine-readable table.
    tsv = args.search_dir / "consensus_topology.tsv"
    with open(tsv, "w") as fh:
        fh.write("sample\treplicon\ttopology\tlength\tunitigs\thairpin\tgfa\n")
        for r in results:
            fh.write(f"{r['sample']}\t{r['replicon']}\t{r['topology']}\t"
                     f"{r['length']}\t{r['unitigs']}\t{int(r['hairpin'])}\t"
                     f"{r['gfa']}\n")
    print(f"\nReport written to : {txt}")
    print(f"Table written to  : {tsv}")


# ── per-cluster mode ──────────────────────────────────────────────────────────

def run_per_cluster(args) -> None:
    qc_dirs = find_qc_dirs(args.search_dir)
    if not qc_dirs:
        print(f"No '*_autocycler/.../qc_pass' directories found under "
              f"'{args.search_dir}'.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(qc_dirs)} qc_pass director"
          f"{'y' if len(qc_dirs) == 1 else 'ies'} "
          f"(gfa: {args.gfa_name}, strict={args.strict}):\n")
    for d in qc_dirs:
        print(f"  {d}")
    print()

    results, skipped, failed_cmd = [], 0, 0
    png_name = Path(args.gfa_name).with_suffix(".png").name

    for qc_path in qc_dirs:
        sample = sample_name(qc_path)
        print("─" * 60)
        print(f"Sample dir : {qc_path}  (sample: {sample})")

        cluster_dirs = sorted(qc_path.glob("cluster_*"))
        if not cluster_dirs:
            print("  (no cluster_* subdirectories found)")
            continue

        for cluster_dir in cluster_dirs:
            cluster = cluster_dir.name
            gfa = cluster_dir / args.gfa_name
            png = cluster_dir / png_name

            if not gfa.exists():
                print(f"  [{cluster}] SKIP – {gfa} not found")
                skipped += 1
                continue

            if not args.no_dotplot:
                if not run_dotplot(gfa, png):
                    failed_cmd += 1
            if args.make_dotplot_only:
                continue

            c = classify_cluster(gfa, strict=args.strict)
            flag = ""
            if c.get("review"):
                flag = "  ⚠ REVIEW"
            elif c["topology"] == "mixed":
                flag = "  *"
            print(f"  [{cluster}] → {c['topology'].upper():8s}  "
                  f"({c['detail']}){flag}")
            results.append(dict(sample=sample, cluster=cluster,
                                png=str(png), gfa=str(gfa), **c))

    if args.make_dotplot_only:
        print(f"\n{'='*60}\nPNG generation complete (no classification).")
        return

    print(f"\n{'='*60}\nSUMMARY\n{'='*60}")
    if not results:
        print("No results to show.")
    else:
        col = max(len(r["sample"]) for r in results)
        header = (f"{'Sample':<{col}}  {'Cluster':<12}  {'Topology':<10}  "
                  f"{'Circular':>10}  Detail")
        print(header)
        print("-" * len(header))
        for r in results:
            frac = f"{r['n_circular']}/{r['n_seqs']}" if r["n_seqs"] else "-"
            print(f"{r['sample']:<{col}}  {r['cluster']:<12}  "
                  f"{r['topology'].upper():<10}  {frac:>10}  {r['detail']}")
        counts = Counter(r["topology"] for r in results)
        print("\nTotals: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
        n_review = sum(1 for r in results if r.get("review"))
        if counts.get("mixed"):
            print("  * mixed = assemblies disagree on topology.")
        if n_review:
            print(f"  ⚠ {n_review} cluster(s) flagged for REVIEW "
                  f"(P-line vote disagrees with graph structure).")

    if skipped:
        print(f"\nSkipped (no GFA): {skipped} cluster(s)")
    if failed_cmd:
        print(f"Failed autocycler dotplot calls: {failed_cmd} "
              f"(classification uses the GFA, not the PNG)")

    if results:
        tsv = args.search_dir / "cluster_topology.tsv"
        with open(tsv, "w") as fh:
            fh.write("sample\tcluster\ttopology\tn_seqs\tn_circular\tn_linear\t"
                     "hairpin\tgraph_topology\treview\tgfa\n")
            for r in results:
                fh.write(f"{r['sample']}\t{r['cluster']}\t{r['topology']}\t"
                         f"{r['n_seqs']}\t{r['n_circular']}\t{r['n_linear']}\t"
                         f"{int(r['hairpin'])}\t{r.get('graph_topology','')}\t"
                         f"{int(bool(r.get('review')))}\t{r['gfa']}\n")
        print(f"\nResults written to: {tsv}")


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify Autocycler clusters/replicons as circular/linear "
                    "from the GFA graph.")
    parser.add_argument("--search-dir", type=Path, default=Path("."),
                        help="Root directory to search from (default: cwd)")
    parser.add_argument("--gfa-name", default="1_untrimmed.gfa",
                        help="Per-cluster GFA filename to classify "
                             "(default: 1_untrimmed.gfa; e.g. 5_final.gfa)")
    parser.add_argument("--strict", action="store_true",
                        help="Also compute graph-structure topology "
                             "(component_is_circular_loop) and flag clusters "
                             "where it disagrees with the P-line vote.")
    parser.add_argument("--consensus", action="store_true",
                        help="Classify autocycler_out/consensus_assembly.gfa "
                             "per replicon (connected component) instead of "
                             "per cluster.")
    parser.add_argument("--no-dotplot", action="store_true",
                        help="Do not run `autocycler dotplot` (skip PNGs).")
    parser.add_argument("--make-dotplot-only", action="store_true",
                        help="Only (re)generate PNGs, do not classify.")
    args = parser.parse_args()

    if args.consensus:
        run_consensus(args)
    else:
        run_per_cluster(args)


if __name__ == "__main__":
    main()

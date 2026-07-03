# linear-plasmid-hairpin-tools

Small, dependency-free Python tools for working out the **topology** of
bacterial replicons from [Autocycler](https://github.com/rrwick/Autocycler)
long-read assemblies — deciding whether each cluster/replicon is **circular** or
**linear**, and detecting the **hairpin / terminal-inverted-repeat ends** that
mark linear plasmids.

These were built while hunting putative linear plasmids in ONT bacterial
assemblies. They read Autocycler and raw-assembler GFA/FASTA files directly and
classify topology from the graph structure and sequence, rather than from
dot-plot pixels (which conflate strand orientation and inverted repeats with
topology).

Three tools:

| Script | What it does |
| --- | --- |
| `autocycler_dotplot_classify.py` | Classify each Autocycler cluster / consensus replicon as **circular / linear / mixed** from the GFA graph. |
| `find_hairpins.py` | Detect **hairpin / terminal-inverted-repeat** ends per assembler (GFA links authoritative; FASTA terminal self-RC heuristic, localised terminal-vs-internal), and optionally extract the flagged contigs. |
| `add_hairpin_edges.py` | Add Autocycler-style **hairpin links** (`L n + n -`) to a raw assembler GFA when a terminal fold-back is present, so the hairpin is visible to Bandage / topology tools. |

---

## Background: how topology is encoded in an Autocycler GFA

Autocycler writes a strand-aware unitig graph. In the GFA:

```
S  <id>  <seq>  <tags...>                     # unitig segment
L  <a> <a_strand> <b> <b_strand>  0M          # link between unitig ends
P  <seq_id> <u1><s1>,<u2><s2>,...  *  LN:i:.. FN:Z:<file> HD:Z:<header>
```

* A **circularising** link is a same-strand self-link — `L n + n +` (and its
  mirror `L n - n -`). Autocycler's `is_isolated_and_circular`.
* A **hairpin** link is an opposite-strand self-link — `L n + n -`
  (3′ fold, `hairpin_end`) or `L n - n +` (5′ fold, `hairpin_start`). It is
  created during `autocycler compress` when an input contig runs into the
  reverse complement of itself, and survives trim → resolve → combine into
  `consensus_assembly.gfa`.
* Each `P` line is one input assembly's contig as an ordered unitig path; a
  contig is **circular** iff its path closes into a loop (a link from the path
  end back to its start).

Raw assemblers (flye, miniasm, raven, …) usually **bake** a fold-back into the
contig sequence instead of adding a self-inverting edge, so their GFAs typically
show 0 hairpin links even when the sequence clearly folds. That difference is
why `add_hairpin_edges.py` exists.

---

## Requirements

* **Python ≥ 3.8** — standard library only, no pip packages.
* Optional: **[Autocycler](https://github.com/rrwick/Autocycler)** on `PATH`, if
  you want `autocycler_dotplot_classify.py` to also render dot-plot PNGs
  (classification itself never needs them).

```bash
git clone https://github.com/<your-user>/linear-plasmid-hairpin-tools.git
cd linear-plasmid-hairpin-tools
chmod +x *.py            # optional
python3 autocycler_dotplot_classify.py --help
```

---

## 1. `autocycler_dotplot_classify.py` — circular vs linear

Classifies replicons from the GFA graph (not from PNG pixels).

**Two input modes**

* **Per-cluster** (default): scans
  `*_autocycler/autocycler_out/clustering/qc_pass/cluster_*/<gfa>` and votes
  across the `P`-line input assemblies (`--gfa-name` selects `1_untrimmed.gfa`
  or `5_final.gfa`).
* **Consensus** (`--consensus`): scans
  `*_autocycler/autocycler_out/consensus_assembly.gfa`, which has no `P` lines,
  and classifies each connected component (replicon) from the graph structure —
  a simple circular loop → `circular`, else `linear`, or `fragmented`.

**Usage**

```bash
# Per-cluster, untrimmed graph, cross-checked against graph structure
python3 autocycler_dotplot_classify.py --search-dir . --gfa-name 1_untrimmed.gfa --strict

# Per-cluster, trimmed/resolved graph (cleaner), no PNGs
python3 autocycler_dotplot_classify.py --search-dir . --gfa-name 5_final.gfa --strict --no-dotplot

# Final combined assembly, per replicon (recommended, most reliable)
python3 autocycler_dotplot_classify.py --search-dir . --consensus
```

**Key options**

* `--search-dir DIR` — root to search (default `.`).
* `--gfa-name NAME` — per-cluster GFA to read (default `1_untrimmed.gfa`).
* `--consensus` — classify `consensus_assembly.gfa` per replicon instead.
* `--strict` — also run Autocycler's `component_is_circular_loop` on the whole
  cluster graph and flag clusters where it disagrees with the `P`-line vote
  (`⚠ REVIEW`).
* `--no-dotplot` / `--make-dotplot-only` — skip or only produce PNGs.

**Outputs**

* Per-cluster: `cluster_topology.tsv`.
* Consensus: `consensus_topology.txt` (the formatted console block) and
  `consensus_topology.tsv` (one row per replicon).

**Recommendation:** `consensus_assembly.gfa` is the most reliable source — each
replicon is a resolved single unitig or clean loop, so topology is unambiguous
and independent of assembler strand disagreement.

---

## 2. `find_hairpins.py` — which assembler produces hairpins

Detects hairpin / terminal-inverted-repeat ends across a directory of assembler
outputs and reports a per-assembler verdict.

**Evidence by file type**

* **GFA (authoritative):** an opposite-strand self-link `L n + n -`. The report
  also prints per-GFA link counts, so `0 hairpins` is confirmed real (not a
  parser miss).
* **FASTA (heuristic, localised):** within a terminal window it finds k-mers
  whose reverse complement also occurs, locates the symmetry axis (fold apex),
  and reports whether the outer arm reaches the contig end. Reaches the end →
  **TERMINAL** hairpin; otherwise **INTERNAL** inverted repeat (e.g. rRNA),
  which is excluded from the verdict.

**Usage**

```bash
# Scan a directory of assembler outputs (canu_*.fasta, flye_*.gfa, ...)
python3 find_hairpins.py --dir .

# Extract the flagged terminal-hairpin contigs for dotplot/alignment
python3 find_hairpins.py --dir . --extract-hairpins

# GFA links only (skip the FASTA heuristic)
python3 find_hairpins.py --dir . --gfa-only
```

**Key options**

* `--dir DIR` — directory to scan (default `.`).
* `-k` (default 31), `--window` (terminal window bp, default 50000).
* `--min-shared` — min supporting reverse-complement k-mers (default 25; raise
  to ~200 to drop small incidental repeats).
* `--edge-tol` — max gap (bp) from the terminus to call a fold TERMINAL
  (default 100).
* `--extract-hairpins` — write flagged terminal-hairpin contigs to
  `hairpin_contigs.fasta` with TIR coordinates in the header, plus
  `hairpin_tir_coords.tsv`.

**Outputs**

* `hairpin_report.txt` (grouped, human-readable) and `hairpin_report.tsv`.
* With `--extract-hairpins`: `hairpin_contigs.fasta` and
  `hairpin_tir_coords.tsv`.

Files are grouped by assembler (prefix before the first `_<number>`), giving a
per-assembler HAIRPIN/clean verdict across replicates.

> **Caveat:** the FASTA method flags any terminal inverted repeat, not only true
> hairpin telomeres. A feature reproduced on the same contig by many assemblers
> is likely real biology; confirm borderline calls with a self dot-plot.

---

## 3. `add_hairpin_edges.py` — annotate a raw GFA with hairpin links

Detects a terminal fold-back on each segment of a raw assembler GFA and appends
the same hairpin link Autocycler uses:

```
3' fold-back  ->  L <seg> + <seg> - 0M
5' fold-back  ->  L <seg> - <seg> + 0M
```

This makes the hairpin visible to Bandage, Autocycler's `topology()`, and
`find_hairpins.py` / the classifier — matching Autocycler's link convention.

**Usage**

```bash
# Annotate one or more GFAs; writes <stem>.hairpins.gfa next to each
python3 add_hairpin_edges.py plassembler_*.gfa --overlap

# Combine several runs into the same report/summary without clobbering
python3 add_hairpin_edges.py more_*.gfa --overlap --append
```

**Key options**

* `-k`, `--window`, `--min-shared`, `--edge-tol` — same detector knobs as
  `find_hairpins.py`.
* `--overlap` — emit the measured arm length as the link overlap (e.g. `4015M`)
  instead of Autocycler's `0M` (better for Bandage; `0M` matches Autocycler).
* `--suffix` (default `.hairpins`) / `--inplace` — output naming.
* `--report PATH` / `--summary PATH` — TSV locations.
* `--append` — append to existing TSVs across runs (header written once).

**Outputs**

* `<stem>.hairpins.gfa` per input (original untouched unless `--inplace`).
* `hairpin_edges.tsv` — one row per added link (segment, strands, overlap, end,
  arm length, gap, support).
* `hairpin_summary.tsv` — one row per GFA (including files with 0 hairpins).

> **Scope / honesty:** this is a topological **annotation**. Autocycler also
> collapses the duplicated fold-back arm into a single unitig so the `0M` link
> reconstructs the fold; here the segment sequence still contains both arms.
> The link is correct for detection/visualisation, not for a length-accurate
> assembly-graph rebuild.

---

## Typical workflow

```bash
# 1. Classify replicons in the final assemblies
python3 autocycler_dotplot_classify.py --search-dir AUTOCYCLER_OUT --consensus

# 2. See which assemblers expose the hairpin on the raw outputs
python3 find_hairpins.py --dir assemblies --extract-hairpins

# 3. Annotate the raw GFAs so the hairpins show up in Bandage / topology tools
python3 add_hairpin_edges.py assemblies/*.gfa --overlap
```

A replicon that comes back `linear` / `fragmented (hairpin)` from tool 1 **and**
shows a reproducible terminal inverted repeat across assemblers in tool 2 is a
strong **linear-plasmid-with-hairpin-ends** candidate.

---

## License

GPL-3.0 — see [LICENSE](LICENSE). Chosen to match Autocycler's own license.

## Acknowledgements

Built around the GFA conventions of
[Autocycler](https://github.com/rrwick/Autocycler) by Ryan Wick. These tools are
independent helpers and are not part of, nor endorsed by, the Autocycler
project.

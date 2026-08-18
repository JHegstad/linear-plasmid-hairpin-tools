#!/usr/bin/env python3
"""
autocycler_clean.py

A standalone Python reimplementation of `autocycler clean`
(https://github.com/rrwick/Autocycler/wiki/Autocycler-clean), for cleaning up
unresolved clusters in Autocycler consensus_assembly.gfa files without needing
the Rust `autocycler` binary installed. The graph model and every algorithm
below (dead-end-safe removal, exclusive-link path merging, circular-loop
closure, tig duplication) were derived by reading Autocycler's own source
(clean.rs, unitig_graph.rs, unitig.rs, graph_simplification.rs, gfa2fasta.rs)
line by line, not guessed from the wiki description alone.

For each input GFA it:
  1. Loads the unitig graph (S/L lines; DP:f: depth tag; CL:Z: colour tag)
     into the same bidirected forward/reverse next/prev link model Autocycler
     uses internally (see the Unitig class below).
  2. [-m/--min-depth] Removes tigs with depth <= the threshold, but ONLY where
     doing so would not create a dead end (a neighbour losing its only
     remaining connection on that side) -- this is exactly Autocycler's own
     `-m` rule, ported line-for-line from `remove_low_depth_unitigs`.
  3. Merges all resulting non-branching ("linear") paths back into single
     tigs, including closing simple circular loops into one self-linked tig
     and preserving hairpin (fold-back) ends -- Autocycler's
     `merge_linear_paths`, reimplemented as a simpler pairwise-confluent
     iteration (see the docstring on merge_linear_paths for why that's
     equivalent to Autocycler's batch version).
  4. [-R/--force-resolve, optional] For any cluster still unresolved after
     step 3, forces it down to a single path/tig: walks outward from the
     highest-depth tig, always taking the highest-depth branch, discarding
     everything else; falls back to duplicating a two-link repeat tig (the
     wiki's manual `-d` move) when a cluster is a terminal-inverted-repeat
     tangle that pruning alone can't resolve. This is NOT dead-end-safe --
     unlike -m, it always fully resolves, even if that means discarding a
     lower-depth alternative that happened to be real biology rather than
     noise. See the "SCOPE / honesty" note on resolve_clusters_by_depth.
  5. Renumbers tigs (by length desc, then sequence, then depth desc, matching
     Autocycler's `renumber_unitigs`).
  6. Writes a cleaned GFA, and optionally a FASTA (like `autocycler gfa2fasta`,
     including its circular=/topology= FASTA header tags).

SCOPE / honesty
---------------
This is an independent reimplementation for use without the Rust toolchain,
verified against Autocycler's algorithms and spot-checked against the wiki's
own worked examples (hairpin-end and terminal-inverted-repeat linear
plasmids), but it is not the upstream binary: exotic graph topologies that
never came up during development may behave subtly differently, particularly
around Autocycler's batched Weak-reference cleanup during merges (this
implementation sidesteps that by fully committing each merge/removal before
considering the next, which is mathematically equivalent for well-formed
graphs but was not exhaustively fuzz-tested against the Rust source). Treat
`-m` output with the same confidence as upstream `autocycler clean -m`;
treat `-R`/--force-resolve output as a strong, automatic starting point that
still deserves a spot-check in Bandage before you trust it on a sample where
the discarded branch had non-trivial depth.

Usage
-----
    Single file:
      ./autocycler_clean.py -i sample_autocycler.gfa -o sample_cleaned.gfa \
          -m 1.0 --fasta

    Batch (all *_autocycler.gfa files in a directory):
      ./autocycler_clean.py --in-dir . --pattern '*_autocycler.gfa' \
          --out-dir cleaned -m 1.0 --fasta

    Force every cluster to fully resolve (adds highest-depth-path pruning +
    terminal-inverted-repeat duplication on top of the safe -m removal):
      ./autocycler_clean.py --in-dir . --out-dir cleaned -m 1.0 -R --fasta
"""

import argparse
import glob
import os
import sys
from typing import Dict, List, Optional, Set, Tuple

Strand = bool  # True = '+', False = '-'
Link = Tuple[int, Strand]

COLOUR = {
    'Consentig': 'steelblue',
    'Anchor': 'forestgreen',
    'Bridge': 'pink',
    'Other': 'orangered',
}
COLOUR_TO_TYPE = {v: k for k, v in COLOUR.items()}

COMPLEMENT = bytes.maketrans(b'ACGTacgtNn', b'TGCAtgcaNn')


def revcomp(seq: bytes) -> bytes:
    return seq.translate(COMPLEMENT)[::-1]


class Unitig:
    # Bidirected-graph link model, ported straight from Autocycler's unitig.rs:
    # a tig has a '+' (forward) end and a '-' (reverse, i.e. reverse-complement)
    # end. forward_next/reverse_next = who this end connects OUT to;
    # forward_prev/reverse_prev = who connects IN to this end. Every physical
    # GFA link is stored on both endpoints (once as a next, once as a prev),
    # and the GFA itself carries both the link and its reverse-complement
    # mirror as separate L lines, so no mirroring needs to happen at parse
    # time (see load_gfa) -- it's already explicit in Autocycler's output.
    __slots__ = ('number', 'seq', 'depth', 'utype',
                 'forward_next', 'forward_prev', 'reverse_next', 'reverse_prev')

    def __init__(self, number: int, seq: bytes, depth: float, utype: str = 'Other'):
        self.number = number
        self.seq = seq
        self.depth = depth
        self.utype = utype
        self.forward_next: List[Link] = []
        self.forward_prev: List[Link] = []
        self.reverse_next: List[Link] = []
        self.reverse_prev: List[Link] = []

    def length(self) -> int:
        return len(self.seq)

    def get_seq(self, strand: Strand) -> bytes:
        return self.seq if strand else revcomp(self.seq)

    def is_anchor(self) -> bool:
        return self.utype == 'Anchor'

    def is_consentig(self) -> bool:
        return self.utype == 'Consentig'

    def is_isolated_and_circular(self) -> bool:
        return (len(self.forward_next) == 1 and len(self.forward_prev) == 1 and
                self.forward_next[0] == (self.number, True) and
                self.forward_prev[0] == (self.number, True))


class Graph:
    def __init__(self):
        self.unitigs: Dict[int, Unitig] = {}
        self.header_extra = ''  # e.g. "\tKM:i:15" if present in input

    def max_unitig_number(self) -> int:
        return max(self.unitigs) if self.unitigs else 0


def quit_with_error(msg: str):
    sys.exit(f"Error: {msg}")


def load_gfa(path: str) -> Graph:
    graph = Graph()
    link_lines = []
    with open(path) as f:
        for raw_line in f:
            line = raw_line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            tag = parts[0]
            if tag == 'H':
                for p in parts[1:]:
                    if p.startswith('KM:i:'):
                        graph.header_extra = '\t' + p
            elif tag == 'S':
                if len(parts) < 3:
                    quit_with_error(f"{path}: segment line does not have enough parts")
                number = int(parts[1])
                seq = parts[2].encode()
                depth = None
                utype = 'Other'
                for p in parts[3:]:
                    if p.startswith('DP:f:'):
                        depth = float(p[5:])
                    elif p.startswith('CL:Z:'):
                        utype = COLOUR_TO_TYPE.get(p[5:], 'Other')
                if depth is None:
                    quit_with_error(
                        f"{path}: could not find a depth tag (e.g. DP:f:10.00) on segment {number}.\n"
                        "Are you sure this is an Autocycler-generated GFA file?")
                graph.unitigs[number] = Unitig(number, seq, depth, utype)
            elif tag == 'L':
                link_lines.append(parts)
            # ignore P lines and anything else

    for parts in link_lines:
        if len(parts) < 6 or parts[5] != '0M':
            quit_with_error(
                f"{path}: non-zero overlap found on a GFA link line.\n"
                "Are you sure this is an Autocycler-generated GFA file?")
        seg1, strand1, seg2, strand2 = int(parts[1]), parts[2] == '+', int(parts[3]), parts[4] == '+'
        if seg1 not in graph.unitigs:
            quit_with_error(f"{path}: link refers to nonexistent unitig {seg1}")
        if seg2 not in graph.unitigs:
            quit_with_error(f"{path}: link refers to nonexistent unitig {seg2}")
        u1, u2 = graph.unitigs[seg1], graph.unitigs[seg2]
        (u1.forward_next if strand1 else u1.reverse_next).append((seg2, strand2))
        (u2.forward_prev if strand2 else u2.reverse_prev).append((seg1, strand1))

    return graph


def delete_dangling_links(graph: Graph):
    valid = graph.unitigs.keys()
    for u in graph.unitigs.values():
        u.forward_next = [l for l in u.forward_next if l[0] in valid]
        u.forward_prev = [l for l in u.forward_prev if l[0] in valid]
        u.reverse_next = [l for l in u.reverse_next if l[0] in valid]
        u.reverse_prev = [l for l in u.reverse_prev if l[0] in valid]


def link_exists(graph: Graph, a_num: int, a_strand: Strand, b_num: int, b_strand: Strand) -> bool:
    a = graph.unitigs.get(a_num)
    if a is None:
        return False
    return (b_num, b_strand) in (a.forward_next if a_strand else a.reverse_next)


# ---------------------------------------------------------------------------
# Step 1: dead-end-safe low-depth removal (equivalent to `autocycler clean -m`)
# ---------------------------------------------------------------------------

def remove_low_depth_unitigs(graph: Graph, min_depth: float) -> int:
    """Removes unitigs with depth <= min_depth, but only where doing so would
    not create a dead end. Returns the number removed."""
    if not graph.unitigs:
        return 0
    original_order = list(graph.unitigs.keys())
    removed = 0
    # Snapshot the order once, then walk it back-to-front. Each deletion is
    # applied immediately (not batched), so a later check in this same pass
    # sees the graph as it stood *after* earlier removals -- exactly
    # Autocycler's own behaviour, and it means removing one low-depth tig can
    # legitimately unlock removal of a neighbour that would otherwise have
    # been a dead end.
    for u_num in reversed(original_order):  # reverse order: keep earlier (usually longer) tigs
        u = graph.unitigs.get(u_num)
        if u is None:
            continue
        if u.depth > min_depth:
            continue

        # Dead-end check: would deleting u strip a neighbour down to zero
        # connections on the side it touches u from? forward_next covers u's
        # '+'-end neighbours, forward_prev its '-'-end neighbours -- together
        # that's every direct connection u has, so checking just these two
        # lists (not reverse_next/reverse_prev, which describe the same
        # physical links from the *other* endpoint's perspective) is enough.
        ok_to_delete = True
        for (nn, nstrand) in u.forward_next:
            if nn == u_num:
                continue
            nb = graph.unitigs[nn]
            plist = nb.forward_prev if nstrand else nb.reverse_prev
            if not any(pn != u_num for pn, _ in plist):
                ok_to_delete = False  # u is nb's only connection on this side -> dead end
                break
        if ok_to_delete:
            for (pn, pstrand) in u.forward_prev:
                if pn == u_num:
                    continue
                pb = graph.unitigs[pn]
                nlist = pb.forward_next if pstrand else pb.reverse_next
                if not any(nn != u_num for nn, _ in nlist):
                    ok_to_delete = False
                    break
        if not ok_to_delete:
            continue

        del graph.unitigs[u_num]
        delete_dangling_links(graph)
        removed += 1
    return removed


# ---------------------------------------------------------------------------
# Step 2: merge linear (non-branching) paths, including circular-loop closure
# ---------------------------------------------------------------------------

def get_exclusive_outputs(graph: Graph, u_num: int) -> List[Link]:
    """u's forward_next neighbours, but ONLY if every one of them has u as its
    sole predecessor on the side it's reached from -- i.e. this is the
    "out-degree(u) == in-degree(neighbour) == 1, mutually" check that defines
    a safely-mergeable, non-branching link. Returns [] (not the mismatched
    list) the moment any neighbour turns out to have another predecessor, so
    callers can just test len(...) == 1 to know a clean single link exists."""
    u = graph.unitigs[u_num]
    outputs = []
    for (nn, nstrand) in u.forward_next:
        nb = graph.unitigs[nn]
        prevs = nb.forward_prev if nstrand else nb.reverse_prev
        if not (len(prevs) == 1 and prevs[0] == (u_num, True)):
            return []
        outputs.append((nn, nstrand))
    if any(nn == u_num for nn, _ in outputs):
        return []  # a link back to itself (hairpin/loop) is never "mergeable onward"
    return outputs


def get_exclusive_inputs(graph: Graph, u_num: int) -> List[Link]:
    """Mirror of get_exclusive_outputs, walking u's forward_prev side instead
    -- used when merge_linear_paths is extending a path along a tig's '-'
    strand, where "next" and "prev" swap roles."""
    u = graph.unitigs[u_num]
    inputs = []
    for (pn, pstrand) in u.forward_prev:
        pb = graph.unitigs[pn]
        nexts = pb.forward_next if pstrand else pb.reverse_next
        if not (len(nexts) == 1 and nexts[0] == (u_num, True)):
            return []
        inputs.append((pn, pstrand))
    if any(pn == u_num for pn, _ in inputs):
        return []
    return inputs


def connected_unitigs(graph: Graph, u_num: int) -> Set[int]:
    u = graph.unitigs[u_num]
    conns = set()
    for lst in (u.forward_next, u.forward_prev, u.reverse_next, u.reverse_prev):
        conns.update(n for n, _ in lst)
    return conns


def connected_components(graph: Graph) -> List[List[int]]:
    visited = set()
    components = []
    for start in graph.unitigs:
        if start in visited:
            continue
        stack = [start]
        comp = []
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            stack.extend(connected_unitigs(graph, cur) - visited)
        comp.sort()
        components.append(comp)
    return components


def component_is_circular_loop(graph: Graph, component: List[int]) -> bool:
    if not component:
        return False
    first = component[0]
    num, strand = first, True
    visited = set()
    while True:
        if num in visited:
            return False
        visited.add(num)
        u = graph.unitigs[num]
        if (len(u.forward_next) != 1 or len(u.forward_prev) != 1 or
                len(u.reverse_next) != 1 or len(u.reverse_prev) != 1):
            return False
        num, strand = (u.forward_next[0] if strand else u.reverse_next[0])
        if num == first:
            break
    return len(visited) == len(component)


def fix_circular_loops(graph: Graph) -> Set[int]:
    fixed_starts = set()
    for component in connected_components(graph):
        if component_is_circular_loop(graph, component):
            fixed_starts.add(component[0])
    return fixed_starts


def weighted_mean_depth(graph: Graph, path: List[Link]) -> float:
    total_len = sum(graph.unitigs[n].length() for n, _ in path)
    if total_len == 0:
        return 0.0
    return sum(graph.unitigs[n].depth * graph.unitigs[n].length() for n, _ in path) / total_len


def get_merge_path_depth(graph: Graph, path: List[Link]) -> float:
    for n, _ in path:
        if graph.unitigs[n].is_anchor():
            return graph.unitigs[n].depth
    return weighted_mean_depth(graph, path)


def merge_path(graph: Graph, path: List[Link], new_num: int):
    merged_seq = b''.join(graph.unitigs[n].get_seq(s) for n, s in path)
    first_num, first_strand = path[0]
    last_num, last_strand = path[-1]
    first_u = graph.unitigs[first_num]
    last_u = graph.unitigs[last_num]

    end_to_start_link = link_exists(graph, last_num, last_strand, first_num, first_strand)
    start_flip_link = link_exists(graph, first_num, not first_strand, first_num, first_strand)
    end_flip_link = link_exists(graph, last_num, last_strand, last_num, not last_strand)

    forward_prev = list(first_u.forward_prev if first_strand else first_u.reverse_prev)
    reverse_next = list(first_u.reverse_next if first_strand else first_u.forward_next)
    forward_next = list(last_u.forward_next if last_strand else last_u.reverse_next)
    reverse_prev = list(last_u.reverse_prev if last_strand else last_u.forward_prev)

    depth = get_merge_path_depth(graph, path)
    utype = 'Consentig' if any(graph.unitigs[n].is_anchor() or graph.unitigs[n].is_consentig()
                                for n, _ in path) else 'Other'

    new_u = Unitig(new_num, merged_seq, depth, utype)
    new_u.forward_next = forward_next
    new_u.forward_prev = forward_prev
    new_u.reverse_next = reverse_next
    new_u.reverse_prev = reverse_prev
    graph.unitigs[new_num] = new_u

    for (nn, nstrand) in list(new_u.forward_next):
        tgt = graph.unitigs[nn]
        (tgt.forward_prev if nstrand else tgt.reverse_prev).append((new_num, True))
    for (nn, nstrand) in list(new_u.forward_prev):
        tgt = graph.unitigs[nn]
        (tgt.forward_next if nstrand else tgt.reverse_next).append((new_num, True))
    for (nn, nstrand) in list(new_u.reverse_next):
        tgt = graph.unitigs[nn]
        (tgt.forward_prev if nstrand else tgt.reverse_prev).append((new_num, False))
    for (nn, nstrand) in list(new_u.reverse_prev):
        tgt = graph.unitigs[nn]
        (tgt.forward_next if nstrand else tgt.reverse_next).append((new_num, False))

    if end_to_start_link:
        new_u.forward_next.append((new_num, True))
        new_u.forward_prev.append((new_num, True))
        new_u.reverse_next.append((new_num, False))
        new_u.reverse_prev.append((new_num, False))
    if start_flip_link:
        new_u.reverse_next.append((new_num, True))
        new_u.forward_prev.append((new_num, False))
    if end_flip_link:
        new_u.forward_next.append((new_num, False))
        new_u.reverse_prev.append((new_num, True))

    for n, _ in path:
        del graph.unitigs[n]


def merge_linear_paths(graph: Graph):
    """Collapse every maximal non-branching path (and simple circular loop)
    into a single tig.

    Autocycler's own version finds all mergeable paths in one pass over the
    pristine graph, then applies every merge_path() call in a batch, relying
    on delete_dangling_links() at the very end to clean up references that
    went stale mid-batch (e.g. two paths that turn out to connect to each
    other). Here we instead merge one adjacent pair at a time and call
    delete_dangling_links() after every single merge, then let the outer
    while-loop pick up wherever a merge just created a new mergeable link
    (including a freshly-closed circular loop -- see the merge_path() self-
    link handling). Path/unipath compaction is confluent: for a graph with no
    ambiguity in what "the next mergeable pair" means, the FINAL merged graph
    is the same regardless of which order you merge pairs in. That lets this
    simpler, fully-committed-per-step approach reach the identical end state
    without needing to reason about Autocycler's specific Weak-reference
    batching order.
    """
    fixed_starts = fix_circular_loops(graph)
    fixed_ends: Set[int] = set()  # always empty: we never merge along tracked input-sequence paths

    def cannot_merge_start(num, strand):
        return (strand and num in fixed_starts) or ((not strand) and num in fixed_ends)

    def cannot_merge_end(num, strand):
        return (strand and num in fixed_ends) or ((not strand) and num in fixed_starts)

    progress = True
    while progress:
        progress = False
        for u_num in list(graph.unitigs.keys()):
            if u_num not in graph.unitigs:
                continue
            for strand in (True, False):
                if cannot_merge_end(u_num, strand):
                    continue
                outputs = get_exclusive_outputs(graph, u_num) if strand else get_exclusive_inputs(graph, u_num)
                if len(outputs) != 1:
                    continue
                out_num, out_strand = outputs[0]
                if not strand:
                    out_strand = not out_strand
                if out_num == u_num:
                    continue
                if cannot_merge_start(out_num, out_strand):
                    continue
                new_num = graph.max_unitig_number() + 1
                merge_path(graph, [(u_num, strand), (out_num, out_strand)], new_num)
                delete_dangling_links(graph)
                progress = True
                break
            if progress:
                break


# ---------------------------------------------------------------------------
# Step 3 (optional, --force-resolve): for any cluster still unresolved after
# the dead-end-safe removal + merge, force it down to a single path by always
# following the highest-depth branch and discarding everything else.
# ---------------------------------------------------------------------------

def strongest_path_in_component(graph: Graph, component: List[int]) -> Set[int]:
    """Starting from the highest-depth tig in this connected component, walk
    outward in both directions, always taking the highest-depth candidate at
    any branch. Returns the set of tig numbers on the resulting path (or
    loop) -- everything else in the component is a discarded side-branch.

    SCOPE: `kept` tracks tig *numbers*, not (number, strand) pairs. That's
    fine for the ordinary case of two different tigs competing at a branch.
    But if a single repeat tig connects the SAME junction from both of its
    strands (a terminal inverted repeat), that tig stays "in kept" the moment
    either of its two strand-orientations is visited, so the walk never
    actually decides between the two parallel links and the component stays
    linked with 2 tigs instead of collapsing to 1. resolve_clusters_by_depth
    detects that stall (removed == 0 with clusters still >1 tig) and falls
    back to duplicate_unitig() for exactly this case -- see the wiki's
    "linear plasmid with terminal inverted repeat" example, which needs the
    manual `-d` step for the same reason.
    """
    comp_set = set(component)
    start = max(component, key=lambda n: (graph.unitigs[n].depth, graph.unitigs[n].length(), n))
    kept: Set[int] = set()

    def walk(num: int, strand: Strand):
        while True:
            kept.add(num)
            u = graph.unitigs[num]
            links = [l for l in (u.forward_next if strand else u.reverse_next) if l[0] in comp_set]
            if not links:
                return
            if len(links) == 1:
                nxt = links[0]
            else:
                nxt = max(links, key=lambda l: (graph.unitigs[l[0]].depth, graph.unitigs[l[0]].length(), l[0]))
            if nxt[0] in kept:
                return  # closed a loop, or met the walk coming from the other direction
            num, strand = nxt

    walk(start, True)   # explore off the '+' end of the highest-depth tig
    walk(start, False)  # explore off the '-' end of the highest-depth tig
    return kept


def resolve_clusters_by_depth(graph: Graph) -> int:
    """For every connected component with more than one tig, keep only the
    single highest-depth path through it and delete the rest. Returns the
    number of tigs removed."""
    removed = 0
    for component in connected_components(graph):
        if len(component) <= 1:
            continue
        kept = strongest_path_in_component(graph, component)
        for n in component:
            if n not in kept:
                del graph.unitigs[n]
                removed += 1
    if removed:
        delete_dangling_links(graph)
    return removed


def _create_link_one_way(graph: Graph, a_num: int, a_strand: Strand, b_num: int, b_strand: Strand):
    """Register a link on both endpoints: a's next list gets b, b's prev list
    gets a. Only records the ONE direction given -- create_link (below) is
    what also adds the reverse-complement mirror, matching how every GFA link
    Autocycler writes out is really two of these calls."""
    a = graph.unitigs[a_num]
    (a.forward_next if a_strand else a.reverse_next).append((b_num, b_strand))
    b = graph.unitigs[b_num]
    (b.forward_prev if b_strand else b.reverse_prev).append((a_num, a_strand))


def create_link(graph: Graph, a_num: int, a_strand: Strand, b_num: int, b_strand: Strand):
    """Add a link a(a_strand) -> b(b_strand) plus its reverse-complement
    mirror b(!b_strand) -> a(!a_strand), i.e. what one physical GFA
    connection actually needs on both ends. The one case that must NOT be
    mirrored a second time is a link that is its own mirror -- a==b with
    opposite strands (a self-palindromic hairpin) -- otherwise it would be
    inserted twice.
    """
    _create_link_one_way(graph, a_num, a_strand, b_num, b_strand)
    if not (a_num == b_num and a_strand != b_strand):
        _create_link_one_way(graph, b_num, not b_strand, a_num, not a_strand)


def duplicate_unitig(graph: Graph, unitig_num: int) -> bool:
    """Splits a tig with exactly two non-self links into two half-depth copies,
    each keeping one of the two links -- for repeats (e.g. a terminal inverted
    repeat) that connect the same junction from both strands, matching what
    the wiki calls `autocycler clean -d`. Returns False if not eligible.

    Ported from Autocycler's duplicate_unitig_by_number: any link the
    original tig had back to ITSELF (a hairpin/loop) is copied onto both new
    copies unchanged; the two links to something else are then handed out
    one each, so each copy keeps exactly one route out.
    """
    u = graph.unitigs.get(unitig_num)
    if u is None:
        return False
    non_self = ([l for l in u.forward_next if l[0] != unitig_num] +
                [l for l in u.reverse_next if l[0] != unitig_num])
    if len(non_self) != 2:
        return False

    a_num = graph.max_unitig_number() + 1
    b_num = a_num + 1
    graph.unitigs[a_num] = Unitig(a_num, u.seq, u.depth / 2.0, u.utype)
    graph.unitigs[b_num] = Unitig(b_num, u.seq, u.depth / 2.0, u.utype)

    orig_forward_next = list(u.forward_next)
    orig_reverse_next = list(u.reverse_next)
    del graph.unitigs[unitig_num]
    delete_dangling_links(graph)

    for (nn, nstrand) in orig_forward_next:
        if nn == unitig_num:
            create_link(graph, a_num, True, a_num, nstrand)
            create_link(graph, b_num, True, b_num, nstrand)
    for (nn, nstrand) in orig_reverse_next:
        if nn == unitig_num:
            create_link(graph, a_num, False, a_num, nstrand)
            create_link(graph, b_num, False, b_num, nstrand)

    non_self_links = ([(True, nn, nstrand) for nn, nstrand in orig_forward_next if nn != unitig_num] +
                       [(False, nn, nstrand) for nn, nstrand in orig_reverse_next if nn != unitig_num])
    (a_strand, a_tgt, a_tgt_strand) = non_self_links[0]
    (b_strand, b_tgt, b_tgt_strand) = non_self_links[1]
    create_link(graph, a_num, a_strand, a_tgt, a_tgt_strand)
    create_link(graph, b_num, b_strand, b_tgt, b_tgt_strand)
    return True


def renumber_unitigs(graph: Graph):
    ordered = sorted(graph.unitigs.values(),
                      key=lambda u: (-u.length(), u.seq, -u.depth))
    old_to_new = {u.number: i + 1 for i, u in enumerate(ordered)}
    new_unitigs: Dict[int, Unitig] = {}
    for u in ordered:
        new_num = old_to_new[u.number]
        nu = Unitig(new_num, u.seq, u.depth, u.utype)
        nu.forward_next = [(old_to_new[n], s) for n, s in u.forward_next]
        nu.forward_prev = [(old_to_new[n], s) for n, s in u.forward_prev]
        nu.reverse_next = [(old_to_new[n], s) for n, s in u.reverse_next]
        nu.reverse_prev = [(old_to_new[n], s) for n, s in u.reverse_prev]
        new_unitigs[new_num] = nu
    graph.unitigs = new_unitigs


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_gfa(graph: Graph, path: str):
    with open(path, 'w') as f:
        f.write(f"H\tVN:Z:1.0{graph.header_extra}\n")
        for u in graph.unitigs.values():
            f.write(f"S\t{u.number}\t{u.seq.decode()}\tDP:f:{u.depth:.2f}\tCL:Z:{COLOUR[u.utype]}\n")
        for u in graph.unitigs.values():
            for (nn, nstrand) in u.forward_next:
                f.write(f"L\t{u.number}\t+\t{nn}\t{'+' if nstrand else '-'}\t0M\n")
            for (nn, nstrand) in u.reverse_next:
                f.write(f"L\t{u.number}\t-\t{nn}\t{'+' if nstrand else '-'}\t0M\n")


def save_fasta(graph: Graph, path: str):
    with open(path, 'w') as f:
        for u in graph.unitigs.values():
            circular = u.is_isolated_and_circular()
            topo = 'circular' if circular else 'linear'
            f.write(f">{u.number} length={u.length()} depth={u.depth:.1f} "
                    f"circular={'true' if circular else 'false'} topology={topo}\n")
            f.write(u.seq.decode() + "\n")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def clean_one(in_gfa: str, out_gfa: str, min_depth: Optional[float] = None,
              out_fasta: Optional[str] = None, force_resolve: bool = False) -> dict:
    graph = load_gfa(in_gfa)
    n_before = len(graph.unitigs)
    n_removed = remove_low_depth_unitigs(graph, min_depth) if min_depth is not None else 0
    merge_linear_paths(graph)
    n_after_safe = len(graph.unitigs)

    n_forced = 0
    n_duplicated = 0
    if force_resolve:
        # Iterate: force the highest-depth path in each remaining cluster, then
        # re-merge (a forced pick can itself unlock a further merge/branch). If
        # pruning stalls with clusters still linked (e.g. a repeat connecting
        # the same junction from both strands, like a terminal inverted
        # repeat), duplicate one eligible repeat tig (the wiki's `-d` move)
        # and keep going.
        for _ in range(50):
            n = resolve_clusters_by_depth(graph)
            if n:
                n_forced += n
                merge_linear_paths(graph)
                continue
            unresolved = [c for c in connected_components(graph) if len(c) > 1]
            if not unresolved:
                break
            duplicated = False
            for component in unresolved:
                for num in component:
                    if duplicate_unitig(graph, num):
                        duplicated = True
                        n_duplicated += 1
                        break
                if duplicated:
                    break
            if not duplicated:
                break  # nothing more can be done automatically
            merge_linear_paths(graph)

    n_after = len(graph.unitigs)
    renumber_unitigs(graph)

    os.makedirs(os.path.dirname(out_gfa) or '.', exist_ok=True)
    save_gfa(graph, out_gfa)
    if out_fasta:
        save_fasta(graph, out_fasta)

    unresolved_components = [c for c in connected_components(graph) if len(c) > 1]
    return {
        'tigs_before': n_before,
        'tigs_removed': n_removed,
        'tigs_after_safe': n_after_safe,
        'tigs_forced_removed': n_forced,
        'tigs_duplicated': n_duplicated,
        'tigs_after': n_after,
        'resolved': len(unresolved_components) == 0,
        'unresolved_clusters': [len(c) for c in unresolved_components],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('-i', '--in-gfa', nargs='+', help='One or more input Autocycler GFA files')
    src.add_argument('--in-dir', help='Directory to batch-process (use with --pattern)')
    ap.add_argument('--pattern', default='*_autocycler.gfa',
                     help='Glob pattern for --in-dir mode (default: *_autocycler.gfa)')
    ap.add_argument('-o', '--out-gfa', help='Output GFA path (single-file mode)')
    ap.add_argument('--out-dir', help='Output directory (batch mode / multi-file -i)')
    ap.add_argument('-m', '--min-depth', type=float, default=None,
                     help='Automatically remove tigs with depth <= this value, provided doing so '
                          'would not create a dead end (same rule as `autocycler clean -m`)')
    ap.add_argument('-R', '--force-resolve', action='store_true',
                     help="After the (optional) safe -m removal and merge, force-resolve any "
                          "cluster that is still unresolved: at every branch, keep only the "
                          "highest-depth tig/path and delete the rest, so every sample ends up "
                          "as a single tig per replicon. Unlike -m, this is NOT dead-end-safe -- "
                          "it always fully resolves, discarding lower-depth alternatives.")
    ap.add_argument('--fasta', action='store_true', help='Also write a cleaned FASTA file')
    ap.add_argument('--suffix', default='_cleaned', help='Suffix inserted before the extension (default: _cleaned)')
    args = ap.parse_args()

    if args.min_depth is None and not args.force_resolve:
        ap.error('specify at least one of -m/--min-depth or -R/--force-resolve')

    if args.in_dir:
        in_files = sorted(glob.glob(os.path.join(args.in_dir, args.pattern)))
        if not in_files:
            sys.exit(f"No files matching {args.pattern!r} found in {args.in_dir}")
        out_dir = args.out_dir or os.path.join(args.in_dir, 'cleaned')
    else:
        in_files = args.in_gfa
        if len(in_files) == 1 and args.out_gfa:
            out_dir = None
        else:
            out_dir = args.out_dir or '.'

    results = []
    for in_gfa in in_files:
        stem = os.path.basename(in_gfa)
        for ext in ('.gfa',):
            if stem.endswith(ext):
                stem = stem[:-len(ext)]
        if stem.endswith('_autocycler'):
            stem = stem[:-len('_autocycler')]

        if out_dir is None:
            out_gfa = args.out_gfa
            out_fasta = os.path.splitext(args.out_gfa)[0] + '.fasta' if args.fasta else None
        else:
            out_gfa = os.path.join(out_dir, f"{stem}{args.suffix}.gfa")
            out_fasta = os.path.join(out_dir, f"{stem}{args.suffix}.fasta") if args.fasta else None

        try:
            stats = clean_one(in_gfa, out_gfa, args.min_depth, out_fasta, args.force_resolve)
        except SystemExit as e:
            print(f"[FAIL] {stem}: {e}", file=sys.stderr)
            continue
        stats['sample'] = stem
        results.append(stats)
        flag = 'OK' if stats['resolved'] else f"UNRESOLVED{stats['unresolved_clusters']}"
        forced_note = f" -{stats['tigs_forced_removed']} forced" if args.force_resolve else ""
        if args.force_resolve and stats['tigs_duplicated']:
            forced_note += f" +{stats['tigs_duplicated']} duplicated(TIR)"
        print(f"[{flag:>18}] {stem:30s} tigs {stats['tigs_before']:3d} -> "
              f"(-{stats['tigs_removed']} low-depth{forced_note}) -> {stats['tigs_after']:3d}")

    n_ok = sum(1 for r in results if r['resolved'])
    print(f"\n{n_ok}/{len(results)} samples fully resolved (single tig per replicon).")
    unresolved = [r['sample'] for r in results if not r['resolved']]
    if unresolved:
        print("Still unresolved, may need manual `-r`/`-d`-style review:")
        for s in unresolved:
            print(f"  - {s}")


if __name__ == '__main__':
    main()

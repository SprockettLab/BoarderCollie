"""Summaries of minimap2 PAF alignments.

aligned fraction: query bases covered by at least one alignment / query length.
identity: summed matching bases / summed alignment block lengths.
embedded: some alignment leaves at least flank_bp of the target sequence on
    one side of it (the match lies inside a longer target sequence).
copies: alignments covering at least min_cover of the query, optionally
    excluding the query's own position in its own genome.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass


@dataclass
class PafRecord:
    query: str
    qlen: int
    qstart: int
    qend: int
    target: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    block: int


def read_paf(path) -> list[PafRecord]:
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        out.append(PafRecord(f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8]),
                             int(f[9]), int(f[10])))
    return out


def union_length(intervals) -> int:
    total, cs, ce = 0, None, None
    for a, b in sorted(intervals):
        if cs is None or a > ce:
            if cs is not None:
                total += ce - cs
            cs, ce = a, b
        else:
            ce = max(ce, b)
    return total + (ce - cs if cs is not None else 0)


@dataclass
class AlignmentSummary:
    aligned_fraction: float
    identity: float
    embedded: bool


def summarize(records: list[PafRecord], lengths: dict[str, int], flank_bp: int) -> dict[str, AlignmentSummary]:
    by_q = defaultdict(list)
    for r in records:
        by_q[r.query].append(r)
    out = {}
    for q, length in lengths.items():
        rs = by_q.get(q, [])
        cov = union_length((r.qstart, r.qend) for r in rs) / length if length else 0.0
        blk = sum(r.block for r in rs)
        ident = sum(r.matches for r in rs) / blk if blk else 0.0
        emb = any(r.tstart >= flank_bp or r.tlen - r.tend >= flank_bp for r in rs)
        out[q] = AlignmentSummary(cov, ident, emb)
    return out


def count_copies(records: list[PafRecord], min_cover: float, own_position: dict[str, tuple[str, int, int]] | None = None):
    """Alignments covering >= min_cover of each query; own_position maps a
    query to (sequence, start, end) of its source (1-based), whose overlapping
    alignments are not counted."""
    counts: dict[str, int] = defaultdict(int)
    for r in records:
        if r.qend - r.qstart < min_cover * r.qlen:
            continue
        if own_position and r.query in own_position:
            seq, s, e = own_position[r.query]
            if r.target == seq and r.tstart < e and r.tend > s - 1:
                continue
        counts[r.query] += 1
    return counts

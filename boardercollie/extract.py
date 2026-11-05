"""Cutting flagged regions and drawing host controls and background windows.

For every candidate region (fcs_report.read_candidates):
  - the region itself;
  - controls_per_segment host controls: windows of the same length at random
    positions in sequences of the same genome that have no candidate region,
    redrawn when more than control_max_n_frac of their bases are N (at most
    200 draws each);
  - background_per_segment windows drawn the same way, used only for
    composition (CpG observed/expected and GC), never aligned.
The random generator is seeded (setting "seed"), so the same inputs give the
same controls.
"""

from __future__ import annotations

import bisect
import random
from pathlib import Path

from .composition import cpg_oe, gc_fraction, percentile
from .fasta import IndexedFasta, write_fasta
from .fcs_report import read_candidates

MAX_DRAWS = 200
SEGMENT_COLUMNS = ("id", "type", "seq_id", "start", "end", "length", "n_frac", "draws",
                   "seq_len", "fcs_action", "fcs_div", "fcs_top_tax", "fcs_top_taxid", "source")
COMPOSITION_COLUMNS = ("id", "gc", "cpg_oe", "background_n", "gc_low", "gc_high", "cpg_oe_low", "cpg_oe_high")


class WindowSampler:
    """Random windows of a given length from a fixed set of sequences."""

    def __init__(self, fa: IndexedFasta, names: list[str], rng: random.Random, max_n: float):
        self.fa, self.rng, self.max_n = fa, rng, max_n
        lengths = fa.lengths()
        self.names = sorted(names, key=lambda n: (-lengths[n], n))
        self.lengths = [lengths[n] for n in self.names]
        self._cum: dict[int, list[int]] = {}

    def draw(self, length: int) -> tuple[str, int, int, str, float, int]:
        # Sequences long enough (names are sorted longest first), each weighted
        # by its number of possible start positions.
        if length not in self._cum:
            k = sum(1 for x in self.lengths if x >= length)
            if k == 0:
                raise ValueError(f"no sequence of at least {length} bp without a flagged region")
            cum, t = [], 0
            for i in range(k):
                t += self.lengths[i] - length + 1
                cum.append(t)
            self._cum[length] = cum
        cum = self._cum[length]
        t = cum[-1]
        for draw in range(1, MAX_DRAWS + 1):
            x = self.rng.randrange(t)
            i = bisect.bisect_right(cum, x)
            start = x - (cum[i - 1] if i else 0) + 1
            seq = self.fa.fetch(self.names[i], start, start + length - 1)
            nf = seq.count("N") / length
            if nf <= self.max_n:
                return self.names[i], start, start + length - 1, seq, nf, draw
        raise ValueError(f"{MAX_DRAWS} draws of {length} bp all had more than {self.max_n:.0%} N")


def extract(genome: Path, action_report: Path, taxonomy_report: Path, settings: dict, fai: Path,
            out_fasta: Path, out_tsv: Path, out_composition: Path) -> int:
    fa = IndexedFasta(genome, fai)
    cands = read_candidates(action_report, taxonomy_report)
    flagged_seqs = {c.seq_id for c in cands}
    unknown = flagged_seqs - set(fa.entries)
    if unknown:
        raise ValueError(f"FCS-GX report names sequences not in the genome: {sorted(unknown)[:5]}")
    rng = random.Random(settings["seed"])
    sampler = WindowSampler(fa, [n for n in fa.entries if n not in flagged_seqs], rng, settings["control_max_n_frac"])
    lo_p, hi_p = settings["conflict_background_percentiles"]

    records, rows, comp = [], [], []
    for i, c in enumerate(cands, 1):
        rid = f"R{i}"
        seq = fa.fetch(c.seq_id, c.start, c.end)
        records.append((rid, seq))
        rows.append([rid, "flagged", c.seq_id, c.start, c.end, c.length, f"{seq.count('N') / len(seq):.4f}", 0,
                     c.seq_len, c.fcs_action, c.fcs_div, c.fcs_top_tax,
                     c.fcs_top_taxid if c.fcs_top_taxid is not None else "NA", c.source])
        for k in range(1, settings["controls_per_segment"] + 1):
            name, s, e, cs, nf, draws = sampler.draw(c.length)
            cid = f"{rid}_C{k}"
            records.append((cid, cs))
            rows.append([cid, "host_control", name, s, e, c.length, f"{nf:.4f}", draws, "", "", "", "", "", ""])
        bg_gc, bg_cpg = [], []
        for _ in range(settings["background_per_segment"]):
            *_, bs, _, _ = sampler.draw(c.length)
            bg_gc.append(gc_fraction(bs))
            bg_cpg.append(cpg_oe(bs))
        comp.append([rid, f"{gc_fraction(seq):.4f}", f"{cpg_oe(seq):.4f}", len(bg_gc),
                     f"{percentile(bg_gc, lo_p):.4f}", f"{percentile(bg_gc, hi_p):.4f}",
                     f"{percentile(bg_cpg, lo_p):.4f}", f"{percentile(bg_cpg, hi_p):.4f}"])
    fa.close()

    write_fasta(records, out_fasta)
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(SEGMENT_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(out_composition, "w") as fh:
        fh.write("\t".join(COMPOSITION_COLUMNS) + "\n")
        for r in comp:
            fh.write("\t".join(map(str, r)) + "\n")
    return len(cands)

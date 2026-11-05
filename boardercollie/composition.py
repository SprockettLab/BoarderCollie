"""Sequence composition: GC fraction and CpG observed/expected.

CpG observed/expected = count(CG) x length / (count(C) x count(G)), over A, C,
G and T only. Vertebrate genomes are depleted of CG (values near 0.2); most
bacteria are not (values near 1).
"""

from __future__ import annotations


def gc_fraction(seq: str) -> float:
    acgt = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0


def cpg_oe(seq: str) -> float:
    c, g = seq.count("C"), seq.count("G")
    acgt = sum(seq.count(b) for b in "ACGT")
    if not c or not g:
        return 0.0
    return seq.count("CG") * acgt / (c * g)


def percentile(values: list[float], p: float) -> float:
    """Percentile p (0 to 100) by linear interpolation between ranks."""
    if not values:
        raise ValueError("no values")
    v = sorted(values)
    k = (len(v) - 1) * p / 100
    lo, hi = int(k), min(int(k) + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)

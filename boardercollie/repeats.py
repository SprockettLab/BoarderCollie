"""Comparing tandem-repeat units (rule R2b, repeat family).

A tandem repeat has no fixed start: the unit ACGTT is the same repeat as GTTAC
(a rotation) and as its reverse complement. Two units are compared by aligning
one, end to end, anywhere inside the other written twice (which contains every
rotation), on both strands, and keeping the best result.
"""

from __future__ import annotations

import edlib

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def unit_identity(a: str, b: str, period_tolerance: float) -> float:
    """Identity of repeat units a and b at the best rotation, on either strand.

    identity = 1 - edit distance / length of the longer unit. Returns 0 when
    either unit is empty or their lengths differ by more than period_tolerance
    (as a fraction of the longer unit).
    """
    a, b = a.upper(), b.upper()
    if not a or not b:
        return 0.0
    longer = max(len(a), len(b))
    if abs(len(a) - len(b)) > period_tolerance * longer:
        return 0.0
    best = 0.0
    for target in (b, reverse_complement(b)):
        r = edlib.align(a, target + target, mode="HW", task="distance")
        best = max(best, 1 - r["editDistance"] / longer)
    return max(best, 0.0)

"""Tandem Repeats Finder (TRF) results in its -ngs format.

Each sequence starts with "@name"; each repeat line has start, end, period,
copy number, consensus size, percent matches, percent indels, score, A, C, G,
T, entropy, consensus unit, repeat sequence (and flanks).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .paf import union_length


@dataclass
class TrfSummary:
    cover: float          # fraction of the sequence in any repeat
    main_cover: float     # fraction covered by the repeat covering the most bases
    main_period: int
    main_unit: str


def read_trf_ngs(path, lengths: dict[str, int]) -> dict[str, TrfSummary]:
    reps = defaultdict(list)
    cur = None
    for line in open(path):
        if line.startswith("@"):
            cur = line[1:].strip()
            continue
        f = line.split()
        if cur is None or len(f) < 14:
            continue
        reps[cur].append((int(f[0]) - 1, int(f[1]), int(f[2]), f[13].upper()))
    out = {}
    for name, length in lengths.items():
        r = reps.get(name, [])
        if not r or not length:
            out[name] = TrfSummary(0.0, 0.0, 0, "")
            continue
        main = max(r, key=lambda x: x[1] - x[0])
        out[name] = TrfSummary(union_length((a, b) for a, b, _, _ in r) / length,
                               (main[1] - main[0]) / length, main[2], main[3])
    return out

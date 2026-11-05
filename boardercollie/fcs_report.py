"""Reading FCS-GX reports into the list of regions BoarderCollie tests.

Two FCS-GX outputs are read, both by column name:
  action report (*.fcs_gx_report.txt): one row per region FCS-GX acts on.
  taxonomy report (*.taxonomy.rpt): one row per sequence, or per piece of a
      sequence split at gaps (seq-id "name~~start..end"), with the best-matching
      taxa and a result such as "contaminant(div)".
Every action row becomes a candidate, given the tax-id of its best match from
the taxonomy report. Every contaminant piece of the taxonomy report that no
action row covers also becomes a candidate (FCS-GX named it but took no action).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

ACTION_COLUMNS = ("seq_id", "start_pos", "end_pos", "seq_len", "action", "div", "agg_cont_cov", "top_tax_name")
TAXONOMY_COLUMNS = ("seq-id", "seq-len", "tax-name-1", "tax-id-1", "div-1", "result")


class ReportError(ValueError):
    """An FCS-GX report is missing a column or has an unreadable row."""


@dataclass
class Candidate:
    seq_id: str
    start: int            # 1-based, inclusive
    end: int              # 1-based, inclusive
    seq_len: int
    fcs_action: str       # EXCLUDE, TRIM, FIX, REVIEW, ... or "none" (taxonomy report only)
    fcs_div: str
    fcs_top_tax: str
    fcs_top_taxid: int | None
    source: str           # "action_report" or "taxonomy_report"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def _header(lines: list[str], need: tuple[str, ...], path: Path) -> dict[str, int]:
    for line in lines:
        if line.startswith("#") and not line.startswith("##"):
            cols = line.lstrip("#").rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(cols)}
            missing = [c for c in need if c not in idx]
            if missing:
                raise ReportError(f"{path}: missing columns {missing}")
            return idx
    raise ReportError(f"{path}: no header line")


@dataclass
class TaxonomyPiece:
    seq_id: str
    start: int
    end: int
    tax_name: str
    tax_id: int | None
    div: str
    result: str


def read_taxonomy_report(path: str | Path) -> list[TaxonomyPiece]:
    path = Path(path)
    lines = path.read_text().splitlines()
    c = _header(lines, TAXONOMY_COLUMNS, path)
    pieces = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        sid = f[c["seq-id"]]
        if "~~" in sid:
            name, rng = sid.split("~~", 1)
            a, b = rng.split("..")
            start, end = int(a), int(b)
        else:
            name, start, end = sid, 1, int(f[c["seq-len"]])
        tid = f[c["tax-id-1"]]
        pieces.append(TaxonomyPiece(name, start, end, f[c["tax-name-1"]],
                                    int(tid) if tid.isdigit() else None, f[c["div-1"]], f[c["result"]]))
    return pieces


def read_candidates(action_report: str | Path, taxonomy_report: str | Path) -> list[Candidate]:
    action_report = Path(action_report)
    lines = action_report.read_text().splitlines()
    c = _header(lines, ACTION_COLUMNS, action_report)
    pieces = read_taxonomy_report(taxonomy_report)
    by_seq: dict[str, list[TaxonomyPiece]] = {}
    for p in pieces:
        by_seq.setdefault(p.seq_id, []).append(p)

    out: list[Candidate] = []
    covered: list[tuple[str, int, int]] = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        seq, s, e = f[c["seq_id"]], int(f[c["start_pos"]]), int(f[c["end_pos"]])
        top = f[c["top_tax_name"]]
        # Tax-id of the best match: the overlapping taxonomy piece naming the same taxon.
        overlap = [p for p in by_seq.get(seq, []) if p.start <= e and p.end >= s]
        named = [p for p in overlap if p.tax_name == top]
        taxid = (named or overlap or [None])[0]
        out.append(Candidate(seq, s, e, int(f[c["seq_len"]]), f[c["action"]], f[c["div"]], top,
                             taxid.tax_id if taxid else None, "action_report"))
        covered.append((seq, s, e))

    for p in pieces:
        if not p.result.startswith("contaminant"):
            continue
        if any(seq == p.seq_id and s <= p.end and e >= p.start for seq, s, e in covered):
            continue
        seq_len = max(q.end for q in by_seq[p.seq_id])
        out.append(Candidate(p.seq_id, p.start, p.end, seq_len, "none", p.div, p.tax_name, p.tax_id,
                             "taxonomy_report"))
    return out

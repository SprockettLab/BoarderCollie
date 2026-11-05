"""Writing a genome's decisions.tsv: one row per flagged region, every test result and the decision."""

from __future__ import annotations

from pathlib import Path

from .rules import GenomeResult

COLUMNS = (
    "region_id", "accession", "seq_id", "start", "end", "length", "source",
    "fcs_action", "fcs_div", "fcs_top_tax", "fcs_top_taxid",
    "alignments", "best_aligned_fraction", "align_class",
    "trf_cover", "main_repeat_cover", "main_unit_length", "lineage_copies", "distant_copies",
    "cpg_oe", "cpg_oe_background_low", "cpg_oe_background_high", "coding_fraction", "prokaryotic_rrna",
    "relative_alignment_test", "rule", "conclusion", "action", "evidence_conflict", "conflict_reasons", "note",
)


def _fmt(v) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def rows(result: GenomeResult) -> list[list[str]]:
    out = []
    for d in result.decisions:
        r, x, c = d.region, d.region.repeats, d.region.composition
        aln = ";".join(f"{a.role}:{a.reference}:{a.aligned_fraction:.3f}:{a.identity:.3f}:{'embedded' if a.embedded else 'not_embedded'}"
                       for a in r.alignments) or "none"
        out.append([_fmt(v) for v in (
            r.region_id, r.accession, r.seq_id, r.start, r.end, r.length, r.source,
            r.fcs_action, r.fcs_div, r.fcs_top_tax, r.fcs_top_taxid,
            aln, d.best_aligned_fraction, d.align_class,
            x.trf_cover if x else None, x.main_repeat_cover if x else None,
            len(x.main_unit) if x and x.main_unit else None,
            x.lineage_copies if x else None, x.distant_copies if x else None,
            c.cpg_oe if c else None, c.cpg_oe_background_low if c else None, c.cpg_oe_background_high if c else None,
            c.coding_fraction if c else None, c.prokaryotic_rrna if c else None,
            d.relative_alignment_test, d.rule, d.conclusion, d.action, d.evidence_conflict, d.conflict_reasons, d.note,
        )])
    return out


def write_decisions(result: GenomeResult, path: str | Path, settings_hash: str) -> None:
    """decisions.tsv with comment lines giving the genome status, control results and settings hash."""
    path = Path(path)
    with path.open("w") as fh:
        fh.write(f"# accession: {result.accession}\n")
        fh.write(f"# status: {result.status}\n")
        fh.write(f"# relative-genome alignment test: {result.relative_alignment_test}\n")
        fh.write(f"# host controls host_like: {result.host_controls_passed} of {result.host_controls_total}\n")
        fh.write(f"# foreign controls foreign_like: {result.foreign_controls_passed} of {result.foreign_controls_total}\n")
        fh.write(f"# settings hash: {settings_hash}\n")
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows(result):
            fh.write("\t".join(row) + "\n")

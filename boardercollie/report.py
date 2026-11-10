"""The contamination report and genome summary of one host (TSV files).

contamination_report.tsv: one row per flagged region of every genome of the host:
where it is, what FCS-GX called it, every test result, and whether it was kept or
removed and why (in words).
genome_summary.tsv: one row per genome: status, counts, removed bases, checksums,
the relative and reference ladder used.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

RULE_MEANING = {
    "2": "kept: FCS-GX's best match is an endogenous retrovirus (a taxon in keep_taxa); not tested",
    "3": "kept: aligns to related host genomes",
    "4": "kept: tandem repeat with copies in the same lineage and none in the most distant genome",
    "5": "kept: tandem repeat whose repeat unit matches a repeat kept by rule 3 or 4",
    "6": "removed: no check kept it",
    "too_short": "FCS-GX's own action: shorter than min_tested_length, too short for the alignment test",
    "fcs_only": "FCS-GX's own action: the relative-genome alignment test was not applied",
}
REGION_COLUMNS = (
    "host", "accession", "organism", "region_id", "seq_id", "start", "end", "length", "source",
    "fcs_action", "fcs_div", "fcs_top_tax", "best_aligned_fraction", "align_class",
    "rule", "reason", "conclusion", "action", "relative_alignment_test",
    "evidence_conflict", "conflict_reasons", "note", "decisions_file",
)
GENOME_COLUMNS = (
    "host", "accession", "organism", "status", "relative_alignment_test", "host_controls",
    "regions_flagged", "regions_kept", "regions_removed", "removed_bp", "original_bp", "reference_bp",
    "original_md5", "reference_md5", "relative", "ladder", "run_key", "fcs_gx",
)


def _screen_source(screen: str | None) -> str:
    """Where a genome's FCS-GX results came from, in words.

    screen: "run:<database build>", "precomputed:<check.json>", or "-"/None (not screened)."""
    if not screen or screen == "-":
        return "not screened"
    kind, _, value = screen.partition(":")
    if kind == "run":
        return f"FCS-GX run by BoarderCollie (database build {value})"
    check = json.loads(Path(value).read_text())
    return (f"reports made earlier: {check['action_report']} and {check['taxonomy_report']} "
            f"(database build {check['database_build'] or 'not stated'}; checked against the genome)")


def _organism(record: dict) -> str:
    rep = record.get("report") or {}
    return (rep.get("organism") or {}).get("organism_name") or f"tax-id {record.get('tax_id', 'NA')}"


def _header(path: Path) -> dict[str, str]:
    out = {}
    for line in open(path):
        if not line.startswith("#"):
            break
        k, _, v = line[1:].strip().partition(": ")
        out[k] = v
    return out


def write_host_report(host: str, genomes: list[dict], regions_out: Path, genomes_out: Path) -> None:
    """genomes: dicts with accession, decisions (path or None), run (run.json path), record
    (download record path), relative, ladder (list), and screen (see _screen_source)."""
    with open(regions_out, "w", newline="") as rf, open(genomes_out, "w", newline="") as gf:
        rw = csv.writer(rf, delimiter="\t", lineterminator="\n")
        gw = csv.writer(gf, delimiter="\t", lineterminator="\n")
        rw.writerow(REGION_COLUMNS)
        gw.writerow(GENOME_COLUMNS)
        for g in genomes:
            record = json.loads(Path(g["record"]).read_text())
            organism = _organism(record)
            run = json.loads(Path(g["run"]).read_text()) if g.get("decisions") else {}
            if not g.get("decisions"):
                gw.writerow([host, g["accession"], organism, "not screened (download_only)", "not applied",
                             "NA", "NA", "NA", "NA", 0, "NA", "NA", record.get("md5", "NA"), record.get("md5", "NA"),
                             "", "", "", _screen_source(g.get("screen"))])
                continue
            head = _header(Path(g["decisions"]))
            rows = list(csv.DictReader((l for l in open(g["decisions"]) if not l.startswith("#")), delimiter="\t"))
            for r in rows:
                rw.writerow([host, g["accession"], organism, r["region_id"], r["seq_id"], r["start"], r["end"],
                             r["length"], r["source"], r["fcs_action"], r["fcs_div"], r["fcs_top_tax"],
                             r["best_aligned_fraction"], r["align_class"], r["rule"],
                             RULE_MEANING.get(r["rule"], r["rule"]), r["conclusion"], r["action"],
                             r["relative_alignment_test"], r["evidence_conflict"], r["conflict_reasons"], r["note"],
                             g["decisions"]])
            gw.writerow([host, g["accession"], organism, head.get("status", "NA"),
                         head.get("relative-genome alignment test", "NA"), head.get("host controls host_like", "NA"),
                         len(rows), sum(r["action"] == "keep" for r in rows), sum(r["action"] == "remove" for r in rows),
                         run.get("removed_bp", "NA"), run.get("original_bp", "NA"), run.get("clean_bp", "NA"),
                         run.get("original_md5", "NA"), run.get("clean_md5", "NA"),
                         g.get("relative") or "", ",".join(g.get("ladder") or []), run.get("run_key", ""),
                         _screen_source(g.get("screen"))])

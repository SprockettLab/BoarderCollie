"""FCS-GX reports made earlier, used in place of running FCS-GX.

A folder given as fcs_gx_reports in the config holds, for each genome or reference
genome, its action report and taxonomy report, named after the genome's key (the
accession, or the stem of a local FASTA's file name):

    <key>.fcs_gx_report.txt  and  <key>.taxonomy.rpt
or as run_gx.py names them:
    <key>.<tax-id>.fcs_gx_report.txt  and  <key>.<tax-id>.taxonomy.rpt

Before a pair is used it is checked against the genome: every sequence the reports
name must be in the genome with the same length, the reports must parse, and the
tax-id the report was run with (when its header states one) must be the genome's.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

from .fcs_report import ACTION_COLUMNS, TAXONOMY_COLUMNS, ReportError, _header, read_candidates

ACTION_SUFFIX = ".fcs_gx_report.txt"
TAXONOMY_SUFFIX = ".taxonomy.rpt"


def find_reports(folder: str | Path, key: str) -> tuple[Path, Path] | None:
    """The (action report, taxonomy report) pair for one genome, or None if the folder has neither.

    Raises ReportError when only one of the two is there, or when several files match.
    """
    folder = Path(folder)
    found = {}
    for suffix in (ACTION_SUFFIX, TAXONOMY_SUFFIX):
        pat = re.compile(re.escape(key) + r"(\.[0-9]+)?" + re.escape(suffix) + "$")
        hits = sorted(p for p in folder.iterdir() if pat.fullmatch(p.name)) if folder.is_dir() else []
        if len(hits) > 1:
            raise ReportError(f"{key}: several FCS-GX reports in {folder}: {[p.name for p in hits]}")
        found[suffix] = hits[0] if hits else None
    action, taxonomy = found[ACTION_SUFFIX], found[TAXONOMY_SUFFIX]
    if action is None and taxonomy is None:
        return None
    if action is None or taxonomy is None:
        have = action or taxonomy
        raise ReportError(f"{key}: {have} has no matching "
                          f"{'taxonomy report' if taxonomy is None else 'action report'} in {folder}")
    return action, taxonomy


def header_info(path: str | Path) -> dict[str, Any]:
    """The run details FCS-GX writes on a report's first line (##[[name, major, minor], {...}])."""
    with open(path) as fh:
        first = fh.readline()
    if not first.startswith("##"):
        return {}
    try:
        data = json.loads(first[2:])
    except json.JSONDecodeError:
        return {}
    meta = data[1] if isinstance(data, list) and len(data) > 1 and isinstance(data[1], dict) else {}
    out: dict[str, Any] = {}

    def walk(d):
        for k, v in d.items():
            if isinstance(v, dict):
                walk(v)
            elif k in ("build-date", "tax-id") and k not in out:
                out[k] = v
    walk(meta)
    return out


def build_date(action: str | Path, taxonomy: str | Path) -> str | None:
    """The FCS-GX database build date the reports were made with, if their headers state it."""
    dates = {str(header_info(p).get("build-date")) for p in (action, taxonomy)} - {"None"}
    if len(dates) > 1:
        raise ReportError(f"{action} and {taxonomy} name different database builds: {sorted(dates)}")
    return dates.pop() if dates else None


def reports_hash(action: str | Path, taxonomy: str | Path) -> str:
    """Short hash of both reports' contents: part of the run key of a genome cleaned with them."""
    h = hashlib.sha256()
    for p in (action, taxonomy):
        h.update(Path(p).read_bytes())
        h.update(b"\0")
    return h.hexdigest()[:12]


def sequence_lengths(fasta: str | Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in lengths:
                    raise ReportError(f"{fasta}: sequence name {name} occurs twice")
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def check_reports(action: str | Path, taxonomy: str | Path, fasta: str | Path, tax_id: int | None) -> dict[str, Any]:
    """Check a precomputed report pair against its genome; return what was checked.

    Raises ReportError on any mismatch: a sequence the genome lacks, a different
    length, a region past the end of its sequence, or a different tax-id.
    """
    action, taxonomy = Path(action), Path(taxonomy)
    lengths = sequence_lengths(fasta)
    problems = []

    lines = action.read_text().splitlines()
    c = _header(lines, ACTION_COLUMNS, action)
    n_action = 0
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        n_action += 1
        seq, end, seq_len = f[c["seq_id"]], int(f[c["end_pos"]]), int(f[c["seq_len"]])
        if seq not in lengths:
            problems.append(f"action report names {seq}, which is not in the genome")
        elif seq_len != lengths[seq] or end > lengths[seq]:
            problems.append(f"action report gives {seq} length {seq_len} (region end {end}); "
                            f"the genome's {seq} is {lengths[seq]} bp")

    lines = taxonomy.read_text().splitlines()
    c = _header(lines, TAXONOMY_COLUMNS, taxonomy)
    n_tax = 0
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        n_tax += 1
        sid = f[c["seq-id"]]
        name, _, rng = sid.partition("~~")
        if name not in lengths:
            problems.append(f"taxonomy report names {name}, which is not in the genome")
            continue
        if rng:
            end = int(rng.split("..")[1])
            if end > lengths[name]:
                problems.append(f"taxonomy report piece {sid} ends past {name}'s {lengths[name]} bp")
        elif int(f[c["seq-len"]]) != lengths[name]:
            problems.append(f"taxonomy report gives {name} length {f[c['seq-len']]}; the genome's is {lengths[name]} bp")

    stated = {header_info(p).get("tax-id") for p in (action, taxonomy)} - {None}
    stated_ids = {int(t) for t in stated if str(t).isdigit()}
    if tax_id is not None and stated_ids and stated_ids != {int(tax_id)}:
        problems.append(f"the reports were run with tax-id {sorted(stated_ids)}; the genome's tax-id is {tax_id}")

    if problems:
        shown = problems[:20]
        more = f" (and {len(problems) - 20} more)" if len(problems) > 20 else ""
        raise ReportError(f"{action.name} does not match {fasta}: " + "; ".join(shown) + more)
    candidates = read_candidates(action, taxonomy)  # parses every row the workflow will read
    return {
        "source": "precomputed",
        "action_report": str(action), "taxonomy_report": str(taxonomy),
        "action_report_sha256_12": _sha(action),
        "taxonomy_report_sha256_12": _sha(taxonomy),
        "reports_hash": reports_hash(action, taxonomy),
        "database_build": build_date(action, taxonomy),
        "tax_id_in_reports": sorted(stated_ids) or "not stated",
        "genome_tax_id": tax_id,
        "action_rows": n_action, "taxonomy_rows": n_tax, "regions_to_test": len(candidates),
        "genome_sequences": len(lengths),
        "checked": "every sequence named is in the genome with the same length; tax-id matches when stated",
    }


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]

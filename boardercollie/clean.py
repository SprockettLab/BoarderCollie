"""Cleaning a genome: the FCS-GX action report filtered to the removals, then checks.

The removals file (decide.py) lists every region to remove. Regions that came
from the action report are passed to FCS-GX's gx clean-genome as their own
action-report rows. Regions FCS-GX named only in its taxonomy report get a row
written here: EXCLUDE when the region is the whole sequence, TRIM when it
touches either end, FIX otherwise. After cleaning, the removed bases must equal
the summed lengths of the removals, and cleaned plus removed bases must equal
the original.
"""

from __future__ import annotations

from pathlib import Path

from .decide import read_tsv
from .fcs_report import ACTION_COLUMNS, _header


def filtered_report(action_report: Path, removals: Path, out: Path) -> int:
    lines = action_report.read_text().splitlines()
    c = _header(lines, ACTION_COLUMNS, action_report)
    rem = read_tsv(removals)
    want = {(r["seq_id"], r["start"], r["end"]) for r in rem if r["source"] == "action_report"}
    found = set()
    with open(out, "w") as fh:
        for line in lines:
            if line.startswith("#"):
                fh.write(line + "\n")
                continue
            f = line.split("\t")
            key = (f[c["seq_id"]], f[c["start_pos"]], f[c["end_pos"]])
            if key in want:
                fh.write(line + "\n")
                found.add(key)
        for r in rem:
            if r["source"] != "taxonomy_report":
                continue
            s, e, n = int(r["start"]), int(r["end"]), int(r["seq_len"])
            action = "EXCLUDE" if (s == 1 and e == n) else "TRIM" if (s == 1 or e == n) else "FIX"
            row = {k: "" for k in ACTION_COLUMNS}
            row.update(seq_id=r["seq_id"], start_pos=s, end_pos=e, seq_len=n, action=action,
                       div="NA", agg_cont_cov="NA", top_tax_name="NA")
            fh.write("\t".join(str(row[k]) for k in ACTION_COLUMNS) + "\n")
    if found != want:
        raise ValueError(f"removal rows not found in {action_report}: {sorted(want - found)[:5]}")
    return len(rem)


def fasta_bases(path: Path) -> tuple[int, int]:
    """(sequences, bases) of a FASTA file."""
    n = b = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
            else:
                b += len(line.rstrip(b"\r\n"))
    return n, b


def check_accounting(original: Path, clean: Path, removed: Path, removals: Path) -> dict:
    want = sum(int(r["length"]) for r in read_tsv(removals))
    on, ob = fasta_bases(original)
    cn, cb = fasta_bases(clean)
    rn, rb = fasta_bases(removed) if removed.exists() else (0, 0)
    if rb != want:
        raise ValueError(f"removed {rb} bp, expected {want} bp (sum of removal lengths)")
    if cb + rb != ob:
        raise ValueError(f"clean {cb} + removed {rb} != original {ob} bp")
    return {"original_sequences": on, "original_bp": ob, "clean_sequences": cn, "clean_bp": cb,
            "removed_sequences": rn, "removed_bp": rb}

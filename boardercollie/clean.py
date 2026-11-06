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

import re
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


def write_removed(original: Path, removals: Path, out: Path, fai: Path) -> int:
    """Every removed region, cut from the original genome (FCS-GX's own
    contamination FASTA holds only whole dropped sequences, not trimmed or
    cut pieces). Returns the bases written."""
    from .fasta import IndexedFasta, write_fasta
    fa = IndexedFasta(original, fai)
    recs = [(f"{r['seq_id']}:{r['start']}-{r['end']} {r['region_id']}",
             fa.fetch(r["seq_id"], int(r["start"]), int(r["end"]))) for r in read_tsv(removals)]
    fa.close()
    write_fasta(recs, out)
    return sum(len(s) for _, s in recs)


GX_LOG = re.compile(r"(\d+) bps dropped; (\d+) bps lowercased; (\d+) bps hardmasked")


def check_accounting(original: Path, clean: Path, removed: Path, removals: Path, gx_log: Path | None = None) -> dict:
    """original bp - clean bp must equal the summed removal lengths, the
    removed FASTA must hold exactly those bases, and FCS-GX must report the
    same number of bases dropped and none lowercased or hardmasked."""
    want = sum(int(r["length"]) for r in read_tsv(removals))
    on, ob = fasta_bases(original)
    cn, cb = fasta_bases(clean)
    rn, rb = fasta_bases(removed) if removed.exists() else (0, 0)
    if ob - cb != want:
        raise ValueError(f"original {ob} - clean {cb} = {ob - cb} bp, expected {want} bp (sum of removal lengths)")
    if rb != want:
        raise ValueError(f"removed FASTA holds {rb} bp, expected {want} bp")
    gx = {}
    if gx_log is not None and gx_log.exists():
        m = GX_LOG.search(gx_log.read_text())
        if not m:
            raise ValueError(f"cannot read the bases dropped, lowercased and hardmasked from {gx_log}")
        dropped, lower, hard = map(int, m.groups())
        if dropped != want or lower or hard:
            raise ValueError(f"FCS-GX reports {dropped} bp dropped, {lower} lowercased, {hard} hardmasked; "
                             f"expected {want} dropped and none lowercased or hardmasked")
        gx = {"gx_bp_dropped": dropped, "gx_bp_lowercased": lower, "gx_bp_hardmasked": hard}
    return {"original_sequences": on, "original_bp": ob, "clean_sequences": cn, "clean_bp": cb,
            "removed_regions": rn, "removed_bp": rb, **gx}

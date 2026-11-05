"""Assembling one genome's evidence from workflow files and applying the rules."""

from __future__ import annotations

import json
from pathlib import Path

from .decisions_io import write_decisions
from .fasta import read_fasta
from .paf import count_copies, read_paf, summarize
from .rules import Alignment, CompositionEvidence, Control, Region, RepeatEvidence, decide_genome
from .trf import read_trf_ngs

REMOVAL_COLUMNS = ("region_id", "seq_id", "start", "end", "length", "seq_len", "fcs_action", "source")


def read_tsv(path) -> list[dict]:
    lines = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    if not lines:
        return []
    head = lines[0]
    return [dict(zip(head, l)) for l in lines[1:]]


def lineage_lookup(taxids, cache: Path, fetch=None):
    """tax-id -> set of ancestors, from the NCBI taxonomy, cached in a JSON file."""
    data = {}
    if cache.exists():
        data = {int(k): set(v) for k, v in json.loads(cache.read_text()).items()}
    missing = sorted(set(taxids) - set(data))
    if missing:
        if fetch is None:
            from .ncbi import taxonomy as fetch
        for k, v in fetch(missing).items():
            data[k] = set(v["parents"])
        cache.write_text(json.dumps({k: sorted(v) for k, v in data.items()}, indent=1))
    return lambda t: data.get(t, set())


def decide_from_files(
    accession: str,
    segments_tsv: Path,
    segments_fasta: Path,
    composition_tsv: Path,
    foreign_fasta: Path,
    alignment_pafs: list[tuple[str, str, Path]],
    trf_ngs: Path,
    own_copies_paf: Path,
    lineage_copies_pafs: list[Path],
    distant_copies_paf: Path,
    settings: dict,
    taxonomy_cache: Path,
    out_decisions: Path,
    out_removals: Path,
    settings_hash: str,
    fetch_taxonomy=None,
):
    segs = read_tsv(segments_tsv)
    flagged = [s for s in segs if s["type"] == "flagged"]
    lengths = {k: len(v) for k, v in read_fasta(segments_fasta).items()}
    foreign = read_fasta(foreign_fasta)
    lengths.update({k: len(v) for k, v in foreign.items()})

    # Alignment summaries per reference.
    summaries = []
    for ref, role, paf in alignment_pafs:
        summaries.append((ref, role, summarize(read_paf(paf), lengths, settings["embedded_flank_bp"])))

    def alignments(qid):
        out = []
        for ref, role, summ in summaries:
            a = summ[qid]
            out.append(Alignment(ref, role, a.aligned_fraction, a.identity, a.embedded if role == "relative" else False))
        return out

    flagged_len = {s["id"]: int(s["length"]) for s in flagged}
    trf = read_trf_ngs(trf_ngs, flagged_len) if flagged else {}
    own_pos = {s["id"]: (s["seq_id"], int(s["start"]), int(s["end"])) for s in flagged}
    cover = settings["repeat_copy_min_cover"]
    own = count_copies(read_paf(own_copies_paf), cover, own_pos)
    lineage = {}
    for p in lineage_copies_pafs:
        for q, n in count_copies(read_paf(p), cover).items():
            lineage[q] = lineage.get(q, 0) + n
    distant = count_copies(read_paf(distant_copies_paf), cover)
    comp = {c["id"]: c for c in read_tsv(composition_tsv)}

    regions = []
    for s in flagged:
        q = s["id"]
        t = trf[q]
        c = comp.get(q)
        taxid = int(s["fcs_top_taxid"]) if s["fcs_top_taxid"] not in ("", "NA") else None
        regions.append(Region(
            q, accession, s["seq_id"], int(s["start"]), int(s["end"]), s["fcs_action"], s["fcs_div"],
            s["fcs_top_tax"], taxid, s["source"], alignments(q),
            RepeatEvidence(t.cover, t.main_cover, t.main_unit, lineage.get(q, 0), distant.get(q, 0)),
            CompositionEvidence(float(c["cpg_oe"]), float(c["cpg_oe_low"]), float(c["cpg_oe_high"])) if c else None,
        ))
    controls = [Control(s["id"], "host", alignments(s["id"])) for s in segs if s["type"] == "host_control"]
    controls += [Control(k, "foreign", alignments(k)) for k in foreign]

    lookup = lineage_lookup([r.fcs_top_taxid for r in regions if r.fcs_top_taxid], taxonomy_cache, fetch_taxonomy)
    result = decide_genome(accession, regions, controls, settings, lookup)
    write_decisions(result, out_decisions, settings_hash)

    seq_len = {s["id"]: s["seq_len"] for s in flagged}
    with open(out_removals, "w") as fh:
        fh.write("\t".join(REMOVAL_COLUMNS) + "\n")
        for d in result.decisions:
            if d.action == "remove":
                r = d.region
                fh.write("\t".join(map(str, (r.region_id, r.seq_id, r.start, r.end, r.length, seq_len[r.region_id],
                                             r.fcs_action, r.source))) + "\n")
    return result

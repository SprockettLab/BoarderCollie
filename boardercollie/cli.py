"""Command-line steps called by the Snakemake workflow (python -m boardercollie.cli <step>)."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
import urllib.request
from pathlib import Path

from . import __version__
from .clean import check_accounting, filtered_report, write_removed
from .config import load_config
from .decide import decide_from_files
from .extract import extract


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_fetch(a):
    from .ncbi import download_genome, genome_report
    out = Path(a.out_fasta)
    rec = {"retrieved": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    if a.local_fasta:
        shutil.copyfile(a.local_fasta, out)
        rec.update(source="local", path=a.local_fasta, tax_id=a.tax_id, md5=md5(out))
    else:
        rec.update(source="ncbi", accession=a.accession, report=genome_report(a.accession))
        rec["md5"] = download_genome(a.accession, out)
        rec["md5_check"] = "matches NCBI md5sum.txt"
        rec["tax_id"] = rec["report"]["organism"]["tax_id"]
    Path(a.out_record).write_text(json.dumps(rec, indent=1))


def cmd_fetch_mitochondrion(a):
    from . import ncbi
    from .mitochondrion import fetch_mitochondrion
    rec = fetch_mitochondrion(a.accession, a.genome, json.loads(Path(a.genome_record).read_text()),
                              Path(a.out_fasta), summary=ncbi.nucleotide_summary, fasta=ncbi.nucleotide_fasta,
                              assembly_mitochondria=ncbi.assembly_mitochondria, ranks_of=ncbi.ranks_of)
    Path(a.out_record).write_text(json.dumps(rec, indent=1))


def cmd_taxid(a):
    print(json.loads(Path(a.record).read_text())["tax_id"])


def cmd_import_fcs_reports(a):
    """Use FCS-GX reports made earlier: check them against the genome, then copy them in."""
    from .fcs_report import ReportError
    from .precomputed import check_reports
    tax_id = json.loads(Path(a.record).read_text()).get("tax_id")
    rec = check_reports(a.action, a.taxonomy, a.fasta, int(tax_id) if tax_id is not None else None)
    if rec["database_build"] and rec["database_build"] != a.build:
        raise ReportError(f"{a.action} names database build {rec['database_build']}; this run uses {a.build}")
    shutil.copyfile(a.action, a.out_report)
    shutil.copyfile(a.taxonomy, a.out_taxonomy)
    rec.update(boardercollie=__version__, genome=a.fasta, genome_md5=md5(Path(a.fasta)), run_database_build=a.build,
               imported=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    Path(a.out_record).write_text(json.dumps(rec, indent=1))


def cmd_foreign(a):
    lines = [l.split("\t") for l in Path(a.tsv).read_text().splitlines() if l and not l.startswith("#")]
    with open(a.out, "w") as fh:
        for acc, *_ in lines:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"
            with urllib.request.urlopen(url, timeout=300) as r:
                text = r.read().decode()
            seq = "".join(text.split("\n")[1:])
            if not text.startswith(">") or not seq:
                raise ValueError(f"no sequence for foreign control {acc}")
            fh.write(f">FOREIGN_{acc}\n{seq}\n")
            time.sleep(0.5)


def cmd_without_flagged(a):
    """Copy a reference genome without every sequence that has any FCS-GX row."""
    from .fcs_report import read_candidates
    drop = {c.seq_id for c in read_candidates(a.report, a.taxonomy)}
    keep, n_drop = True, 0
    with open(a.fasta) as src, open(a.out, "w") as dst:
        for line in src:
            if line.startswith(">"):
                keep = line[1:].split()[0] not in drop
                n_drop += not keep
            if keep:
                dst.write(line)
    Path(a.log).write_text(f"{n_drop} sequences with an FCS-GX row removed: {sorted(drop)}\n")


def cmd_extract(a):
    cfg = load_config(a.config)
    n = extract(Path(a.genome), Path(a.report), Path(a.taxonomy), cfg.settings, Path(a.fai),
                Path(a.out_fasta), Path(a.out_tsv), Path(a.out_composition))
    print(f"{n} candidate regions")


def lineage_members(own_record: Path, candidates: list[str], rank: str, fetch=None):
    """Of candidate "accession=record=paf" entries, the PAFs of assemblies in the same taxon at rank."""
    if fetch is None:
        from .ncbi import taxonomy as fetch
    if not candidates:
        return []
    own = json.loads(own_record.read_text())["tax_id"]
    entries = [c.split("=", 2) for c in candidates]
    taxids = {acc: json.loads(Path(rec).read_text())["tax_id"] for acc, rec, _ in entries}
    tax = fetch([own, *taxids.values()])
    own_rank = tax[own]["ranks"].get(rank)
    if own_rank is None:
        raise ValueError(f"tax-id {own} has no {rank} in the NCBI taxonomy")
    keep = [(acc, Path(paf)) for acc, _, paf in entries if tax[taxids[acc]]["ranks"].get(rank) == own_rank]
    return keep


def cmd_decide(a):
    cfg = load_config(a.config)
    pafs = [tuple(x.split(":", 2)) for x in a.alignment]  # reference:role:path
    pafs = [(r, role, Path(p)) for r, role, p in pafs]
    lineage = lineage_members(Path(a.own_record), a.lineage_copies, cfg.settings["repeat_lineage_rank"])
    Path(a.out_lineage).write_text("".join(f"{acc}\t{p}\n" for acc, p in lineage))
    res = decide_from_files(a.accession, Path(a.segments), Path(a.segments_fasta), Path(a.composition),
                            Path(a.foreign), pafs, Path(a.trf) if a.trf else None,
                            Path(a.own_copies) if a.own_copies else None,
                            [p for _, p in lineage], Path(a.distant_copies) if a.distant_copies else None,
                            cfg.settings, Path(a.taxonomy_cache), Path(a.out_decisions), Path(a.out_removals),
                            a.settings_hash, test_unavailable=a.test_unavailable)
    print(f"{a.accession}: {res.status}; relative-genome alignment test {res.relative_alignment_test}")


def cmd_filter_report(a):
    print(filtered_report(Path(a.report), Path(a.removals), Path(a.out)), "removals")


def cmd_check_clean(a):
    write_removed(Path(a.original), Path(a.removals), Path(a.removed), Path(a.fai))
    acc = check_accounting(Path(a.original), Path(a.clean), Path(a.removed), Path(a.removals),
                           Path(a.gx_log) if a.gx_log else None)
    rec = {"boardercollie": __version__, "accession": a.accession, "run_key": a.run_key,
           "settings_hash": a.settings_hash, "status_file": a.decisions,
           "original_md5": md5(Path(a.original)), "clean_md5": md5(Path(a.clean)),
           "removed_md5": md5(Path(a.removed)) if Path(a.removed).exists() else None, **acc,
           "finished": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    status = [l for l in open(a.decisions) if l.startswith("# ")]
    rec["decisions_header"] = [l[2:].strip() for l in status]
    Path(a.out).write_text(json.dumps(rec, indent=1))


def cmd_combine(a):
    names = set()
    # member_of: "-" for a genome; for a genome's mitochondrion, that genome's accession.
    member_of = getattr(a, "member_of", None) or ["-"] * len(a.accessions)
    if not len(a.accessions) == len(a.fastas) == len(a.records) == len(member_of):
        raise ValueError("combine: --accessions, --fastas, --records and --member-of differ in length")
    with open(a.out, "w") as dst, open(a.members_out, "w") as mem:
        mem.write("accession\tclean_genome\tclean_md5\trun_record\tmember_of\n")
        for acc, fasta, record, of in zip(a.accessions, a.fastas, a.records, member_of):
            with open(fasta) as src:
                for line in src:
                    if line.startswith(">"):
                        n = line[1:].split()[0]
                        if n in names:
                            raise ValueError(f"sequence name {n} occurs in two members")
                        names.add(n)
                    dst.write(line)
            mem.write(f"{acc}\t{fasta}\t{md5(Path(fasta))}\t{record}\t{of}\n")


def cmd_report(a):
    from .report import write_host_report
    genomes = []
    # accession|decisions|run|record|relative|ladder (comma-separated)|screen; "-" for none
    for spec in a.genome:
        acc, dec, run, rec, rel, lad, screen = spec.split("|")
        genomes.append({"accession": acc, "decisions": None if dec == "-" else dec, "run": run, "record": rec,
                        "relative": None if rel == "-" else rel, "ladder": [] if lad == "-" else lad.split(","),
                        "screen": None if screen == "-" else screen})
    write_host_report(a.host, genomes, Path(a.out_regions), Path(a.out_genomes))


def main(argv=None):
    p = argparse.ArgumentParser(prog="boardercollie")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("fetch"); s.set_defaults(f=cmd_fetch)
    s.add_argument("--accession"); s.add_argument("--local-fasta"); s.add_argument("--tax-id", type=int)
    s.add_argument("--out-fasta", required=True); s.add_argument("--out-record", required=True)

    s = sub.add_parser("fetch-mitochondrion"); s.set_defaults(f=cmd_fetch_mitochondrion)
    for k in ("accession", "genome", "genome-record", "out-fasta", "out-record"):
        s.add_argument(f"--{k}", required=True)

    s = sub.add_parser("taxid"); s.set_defaults(f=cmd_taxid); s.add_argument("record")

    s = sub.add_parser("import-fcs-reports", help="check FCS-GX reports made earlier against the genome and copy them")
    s.set_defaults(f=cmd_import_fcs_reports)
    for k in ("fasta", "record", "action", "taxonomy", "build", "out-report", "out-taxonomy", "out-record"):
        s.add_argument(f"--{k}", required=True)

    s = sub.add_parser("foreign-controls"); s.set_defaults(f=cmd_foreign)
    s.add_argument("--tsv", required=True); s.add_argument("--out", required=True)

    s = sub.add_parser("without-flagged"); s.set_defaults(f=cmd_without_flagged)
    for k in ("fasta", "report", "taxonomy", "out", "log"):
        s.add_argument(f"--{k}", required=True)

    s = sub.add_parser("extract"); s.set_defaults(f=cmd_extract)
    for k in ("config", "genome", "report", "taxonomy", "fai", "out-fasta", "out-tsv", "out-composition"):
        s.add_argument(f"--{k}", required=True)

    s = sub.add_parser("decide"); s.set_defaults(f=cmd_decide)
    for k in ("config", "accession", "segments", "segments-fasta", "composition", "foreign",
              "own-record", "taxonomy-cache", "out-decisions", "out-removals",
              "out-lineage", "settings-hash"):
        s.add_argument(f"--{k}", required=True)
    s.add_argument("--trf", default=None, help="TRF output; absent when the alignment checks are off")
    s.add_argument("--own-copies", default=None, help="copies in the genome itself; absent when the alignment checks are off")
    s.add_argument("--distant-copies", default=None, help="copies in the most distant comparison genome; absent when there is none")
    s.add_argument("--test-unavailable", choices=("no_relative", "no_ladder"), default=None,
                   help="why the relative-genome alignment test cannot run; FCS-GX's own actions are applied")
    s.add_argument("--alignment", action="append", default=[], help="reference:role:paf")
    s.add_argument("--lineage-copies", nargs="*", default=[], help="accession=record=paf")

    s = sub.add_parser("filter-report"); s.set_defaults(f=cmd_filter_report)
    for k in ("report", "removals", "out"):
        s.add_argument(f"--{k}", required=True)

    s = sub.add_parser("check-clean"); s.set_defaults(f=cmd_check_clean)
    for k in ("accession", "original", "clean", "removed", "removals", "decisions", "run-key", "settings-hash", "out", "fai"):
        s.add_argument(f"--{k}", required=True)
    s.add_argument("--gx-log", default=None, help="gx clean-genome log; absent when nothing was removed")

    s = sub.add_parser("combine"); s.set_defaults(f=cmd_combine)
    s.add_argument("--accessions", nargs="+", required=True); s.add_argument("--fastas", nargs="+", required=True)
    s.add_argument("--records", nargs="+", required=True)
    s.add_argument("--out", required=True); s.add_argument("--members-out", required=True)
    s.add_argument("--member-of", nargs="+", help="per member: - for a genome, else the genome it belongs to")

    s = sub.add_parser("report"); s.set_defaults(f=cmd_report)
    s.add_argument("--host", required=True)
    s.add_argument("--genome", action="append", required=True,
                   help="accession|decisions|run|record|relative|ladder|screen")
    s.add_argument("--out-regions", required=True); s.add_argument("--out-genomes", required=True)

    a = p.parse_args(argv)
    a.f(a)


if __name__ == "__main__":
    sys.exit(main())

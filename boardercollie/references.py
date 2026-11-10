"""Choosing the closest relative of a genome from the NCBI taxonomy.

The relative is the assembly that shares the deepest common ancestor with the
screened genome: the search starts in the genome's own species, then its genus,
family, order and class, and stops at the first of these that holds any
eligible assembly.

Never eligible:
  - the screened assembly itself, and its GenBank/RefSeq twin (paired accession);
  - any assembly from the same BioSample (the same animal), whose
    contamination would be the same as the screened genome's.

Among the eligible assemblies at that level, the order of preference is:
  1. a different BioProject from the screened genome (two assemblies from one
     sequencing project could share its laboratory contamination);
  2. NCBI's "reference genome" designation;
  3. a more complete assembly level (complete, chromosome, scaffold, contig);
  4. a longer contig N50;
  5. RefSeq (GCF_) over GenBank (GCA_) for the same assembly, then accession,
     so the choice is the same every time for the same NCBI records.

Every candidate considered, and why the chosen one won, is written out so the
choice can be checked.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

SEARCH_RANKS = ("species", "genus", "family", "order", "class")
LEVEL_ORDER = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
AUTOMATIC = {"", "none", "na", "auto", "automatic"}


def is_automatic(value) -> bool:
    """True when a config value asks for automatic choice: missing, blank, none, NA or auto."""
    return value is None or (isinstance(value, str) and value.strip().lower() in AUTOMATIC)


@dataclass
class Assembly:
    accession: str
    paired_accession: str | None
    tax_id: int
    organism: str
    bioproject: str | None
    biosample: str | None
    level: str | None
    reference_genome: bool
    contig_n50: int
    scaffold_n50: int = 0

    @classmethod
    def from_report(cls, r: dict) -> "Assembly":
        a = r.get("assembly_info") or {}
        s = r.get("assembly_stats") or {}
        o = r.get("organism") or {}
        return cls(
            accession=r["accession"],
            paired_accession=r.get("paired_accession"),
            tax_id=int(o.get("tax_id", 0)),
            organism=o.get("organism_name", ""),
            bioproject=a.get("bioproject_accession"),
            biosample=(a.get("biosample") or {}).get("accession"),
            level=a.get("assembly_level"),
            reference_genome=a.get("refseq_category") == "reference genome",
            contig_n50=int(s.get("contig_n50") or 0),
            scaffold_n50=int(s.get("scaffold_n50") or 0),
        )


def _same_assembly(a: Assembly, b: Assembly) -> bool:
    ids_a = {a.accession, a.paired_accession} - {None}
    ids_b = {b.accession, b.paired_accession} - {None}
    return bool(ids_a & ids_b)


def _dedupe(assemblies: list[Assembly]) -> list[Assembly]:
    """One record per assembly: RefSeq (GCF_) kept over its GenBank (GCA_) twin."""
    out: list[Assembly] = []
    for a in sorted(assemblies, key=lambda x: (not x.accession.startswith("GCF_"), x.accession)):
        if not any(_same_assembly(a, b) for b in out):
            out.append(a)
    return out


def preference_key(c: Assembly, screened: Assembly):
    return (
        c.bioproject is not None and c.bioproject == screened.bioproject,
        not c.reference_genome,
        LEVEL_ORDER.get(c.level or "", 9),
        -c.contig_n50,
        not c.accession.startswith("GCF_"),
        c.accession,
    )


def eligible(c: Assembly, screened: Assembly, min_scaffold_n50: int = 0) -> tuple[bool, str]:
    if _same_assembly(c, screened):
        return False, "the screened assembly or its GenBank/RefSeq twin"
    if screened.biosample and c.biosample == screened.biosample:
        return False, "same BioSample (same animal) as the screened genome"
    if c.scaffold_n50 < min_scaffold_n50:
        return False, f"scaffold N50 {c.scaffold_n50} below {min_scaffold_n50}"
    return True, ""


def choose_relative(
    screened: Assembly,
    ranks: dict[str, int],
    genomes_under: Callable[[int], list[dict]],
    min_scaffold_n50: int = 0,
) -> dict:
    """The closest relative of `screened`, with every candidate considered.

    ranks: the screened genome's NCBI classification, rank -> tax-id.
    genomes_under: tax-id -> NCBI assembly reports under it (ncbi.genomes_under).
    Returns {"relative": accession or None, "rank": ..., "taxon": ..., "reason": ..., "candidates": [...]}.
    """
    seen: set[str] = set()
    rows: list[dict] = []
    for rank in SEARCH_RANKS:
        taxon = ranks.get(rank)
        if taxon is None:
            continue
        found = _dedupe([Assembly.from_report(r) for r in genomes_under(taxon)])
        ok = []
        for c in found:
            if c.accession in seen:
                continue
            seen.add(c.accession)
            good, why = eligible(c, screened, min_scaffold_n50)
            rows.append({**asdict(c), "rank": rank, "eligible": good, "excluded_because": why})
            if good:
                ok.append(c)
        if ok:
            ok.sort(key=lambda c: preference_key(c, screened))
            best = ok[0]
            same_project = best.bioproject is not None and best.bioproject == screened.bioproject
            reason = (f"deepest shared rank with any eligible assembly: {rank} (tax-id {taxon}); "
                      f"{len(ok)} eligible there; chosen by: "
                      + ("same BioProject (no other eligible assembly at this rank), " if same_project
                         else "different BioProject, ")
                      + ("NCBI reference genome, " if best.reference_genome else "")
                      + f"assembly level {best.level}, contig N50 {best.contig_n50}")
            for row in rows:
                row["chosen"] = row["accession"] == best.accession
            return {"relative": best.accession, "rank": rank, "taxon": taxon, "reason": reason, "candidates": rows}
    return {"relative": None, "rank": None, "taxon": None,
            "reason": f"no eligible assembly within the genome's {SEARCH_RANKS[-1]}", "candidates": rows}


LADDER_RANKS = ("genus", "family", "order", "class")
RANK_BELOW = {"genus": "species", "family": "genus", "order": "family", "class": "order"}


def choose_ladder(
    host: list[tuple[Assembly, dict[str, int]]],
    exclude: set[str],
    ladder_ranks: list[str],
    reference_genomes_under: Callable[[int], list[dict]],
    lineage_of: Callable[[list[int]], dict[int, dict[str, int]]],
    min_scaffold_n50: int = 0,
) -> dict:
    """The comparison genomes of a host: one per rank in ladder_ranks, nearest first.

    The genome for rank R lies in the host's taxon at R but outside the host's
    taxa at the rank below (for family: same family, other genus), so each one is
    farther away than the one before. Only NCBI reference genomes are considered.
    The last one is the distant reference of rule 4.

    host: (assembly, rank -> tax-id) for each genome of the host.
    exclude: accessions never used (the host's genomes and their relatives).
    reference_genomes_under: tax-id -> NCBI reference genome reports under it.
    lineage_of: tax-ids -> {tax-id: rank -> tax-id}.
    """
    projects = {a.bioproject for a, _ in host if a.bioproject}
    chosen: list[str] = []
    rungs = []
    for rank in ladder_ranks:
        lower = RANK_BELOW[rank]
        taxa = {r.get(rank) for _, r in host}
        if None in taxa:
            rungs.append({"rank": rank, "chosen": None, "reason": f"a host genome has no {rank} in the NCBI taxonomy",
                          "candidates": []})
            continue
        if len(taxa) > 1:
            raise ValueError(f"the host's genomes lie in different taxa at rank {rank} ({sorted(taxa)}); "
                             "name the comparison genomes (ladder) in the config")
        taxon = taxa.pop()
        inner = {r.get(lower) for _, r in host} - {None}
        found = _dedupe([Assembly.from_report(x) for x in reference_genomes_under(taxon)])
        lin = lineage_of(sorted({c.tax_id for c in found})) if found else {}
        rows, ok = [], []
        for c in found:
            why = ""
            if any(_same_assembly(c, Assembly(e, None, 0, "", None, None, None, False, 0)) for e in exclude):
                why = "one of the host's genomes or their relatives"
            elif c.accession in chosen:
                why = "already chosen at a nearer rank"
            elif (lin.get(c.tax_id) or {}).get(lower) in inner:
                why = f"same {lower} as the host"
            elif c.scaffold_n50 < min_scaffold_n50:
                why = f"scaffold N50 {c.scaffold_n50} below {min_scaffold_n50}"
            rows.append({**asdict(c), "eligible": not why, "excluded_because": why})
            if not why:
                ok.append(c)
        if not ok:
            rungs.append({"rank": rank, "taxon": taxon, "chosen": None,
                          "reason": f"no eligible NCBI reference genome in the host's {rank} outside its {lower}",
                          "candidates": rows})
            continue
        ok.sort(key=lambda c: (c.bioproject in projects, *preference_key(c, Assembly("", None, 0, "", None, None, None, False, 0))[1:]))
        best = ok[0]
        chosen.append(best.accession)
        for row in rows:
            row["chosen"] = row["accession"] == best.accession
        rungs.append({"rank": rank, "taxon": taxon, "chosen": best.accession,
                      "reason": f"same {rank}, other {lower}; {len(ok)} eligible; chosen: {best.organism}, "
                                + ("NCBI reference genome, " if best.reference_genome else "")
                                + f"assembly level {best.level}, contig N50 {best.contig_n50}",
                      "candidates": rows})
    return {"ladder": chosen, "rungs": rungs}


def _local_assembly(key: str, tax_id: int) -> Assembly:
    return Assembly(key, None, tax_id, "", None, None, None, False, 0)


def resolve(cfg, library_dir: Path, ncbi=None) -> list[str]:
    """Fill in every relative and ladder the config leaves blank; returns one line per choice.

    Choices are saved in the library (genomes/<genome>/relative_choice.json and
    hosts/<host>/ladder_choice.json) with every candidate considered, and reused
    while the genomes and settings they depend on are unchanged. NCBI is asked only
    when a choice has to be made.
    """
    if ncbi is None:
        from . import ncbi
    s = cfg.settings
    min_n50 = s["relative_min_scaffold_n50"]
    lines: list[str] = []
    info: dict[str, tuple[Assembly, dict[str, int]]] = {}

    def assembly(g):
        if g.key not in info:
            if g.fasta:
                a = _local_assembly(g.key, g.tax_id)
            else:
                a = Assembly.from_report(ncbi.genome_report(g.accession))
            info[g.key] = (a, ncbi.ranks_of([a.tax_id])[a.tax_id])
        return info[g.key]

    for h in cfg.hosts:
        for g in h.genomes:
            if g.relative is not None:
                continue
            inputs = {"genome": g.key, "tax_id": g.tax_id, "relative_min_scaffold_n50": min_n50}
            ch = cached_choice(library_dir / "genomes" / g.key / "relative_choice.json",
                               lambda g=g: choose_relative(*assembly(g), ncbi.genomes_under, min_n50), inputs)
            g.relative = ch["relative"]
            lines.append(f"{h.name} {g.key}: relative {g.relative or 'none found'} ({ch['reason']})")
        if not h.ladder:
            relatives = sorted({g.relative for g in h.genomes} - {None})
            inputs = {"genomes": sorted(g.key for g in h.genomes), "relatives": relatives,
                      "ladder_ranks": s["ladder_ranks"], "relative_min_scaffold_n50": min_n50}
            ch = cached_choice(
                library_dir / "hosts" / h.name / "ladder_choice.json",
                lambda h=h, relatives=relatives: choose_ladder(
                    [assembly(g) for g in h.genomes], {g.key for g in h.genomes} | set(relatives), s["ladder_ranks"],
                    lambda t: ncbi.genomes_under(t, reference_only=True), ncbi.ranks_of, min_n50),
                inputs)
            h.ladder = list(ch["ladder"])
            lines.append(f"{h.name}: comparison genomes {', '.join(h.ladder) or 'none found'} ("
                         + "; ".join(f"{r['rank']}: {r['chosen'] or r['reason']}" for r in ch["rungs"]) + ")")
    return lines


def cached_choice(path: Path, compute: Callable[[], dict], inputs: dict | None = None) -> dict:
    """A choice made once and reused: read from `path` if it exists, otherwise computed and written there.

    NCBI adds assemblies over time, so a fresh search could choose differently; reusing
    the recorded choice keeps a genome's results stable. The choice is made again only
    when `inputs` (the genomes and settings it depends on) differ from those recorded,
    or when the file is deleted.
    """
    inputs = inputs or {}
    if path.exists():
        old = json.loads(path.read_text())
        if old.get("inputs", {}) == inputs:
            return old
    choice = compute()
    choice["inputs"] = inputs
    choice["chosen_on"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(json.dumps(choice, indent=1))
    tmp.rename(path)
    return choice

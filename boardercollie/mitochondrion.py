"""A host genome's mitochondrion, fetched from NCBI Nucleotide and checked.

Most nuclear assemblies are deposited without their mitochondrion: labs usually
submit it as a separate Nucleotide record, or not at all. NCBI's RefSeq copies
(GCF_) carry the species' RefSeq mitochondrial genome (an NC_ record) added after
assembly. A host genome may name such a record (config key "mitochondrion"); it is
added to the host reference right after the genome, under these rules, each checked
here before anything is written:

  1. the record is a complete mitochondrial genome (NCBI: genome "mitochondrion",
     completeness "complete"), and its sequence has the length NCBI states;
  2. it is from the genome's own species (same NCBI species tax-id);
  3. the assembly does not already hold a mitochondrion (an NCBI assembly is asked
     for its sequences; a local FASTA cannot be, which the record says).

The record is not screened with FCS-GX: it is one 15-18 kb curated sequence of the
host's own species, the one NCBI itself adds to RefSeq assemblies.
"""

from __future__ import annotations

import hashlib
import time
from pathlib import Path
from typing import Callable


class MitochondrionError(ValueError):
    """The named record is not a usable mitochondrion for this genome."""


def species_of(taxid: int, ranks: dict[int, dict[str, int]]) -> int:
    sp = (ranks.get(int(taxid)) or {}).get("species")
    if sp is None:
        raise MitochondrionError(f"NCBI taxonomy names no species for tax-id {taxid}")
    return int(sp)


def fetch_mitochondrion(accession: str, genome: str, genome_record: dict, out_fasta: Path, *,
                        summary: Callable[[str], dict], fasta: Callable[[str], str],
                        assembly_mitochondria: Callable[[str], list[str]],
                        ranks_of: Callable[[list[int]], dict[int, dict[str, int]]]) -> dict:
    """Check the record against the genome, write it to out_fasta, return what was checked.

    The NCBI lookups are passed in, so the checks can be tested without the network.
    """
    s = summary(accession)
    if s.get("genome") != "mitochondrion" or s.get("completeness") != "complete":
        raise MitochondrionError(f"{accession} is not a complete mitochondrial genome "
                                 f"(NCBI: genome {s.get('genome')!r}, completeness {s.get('completeness')!r}): "
                                 f"{s.get('title')}")
    mt_taxid, g_taxid = int(s["taxid"]), int(genome_record["tax_id"])
    ranks = ranks_of([mt_taxid, g_taxid])
    if species_of(mt_taxid, ranks) != species_of(g_taxid, ranks):
        raise MitochondrionError(f"{accession} ({s.get('organism')}, tax-id {mt_taxid}) is not from the species of "
                                 f"{genome} (tax-id {g_taxid}); a host's mitochondrion must be its own species'")
    if genome_record.get("source") == "ncbi":
        held = assembly_mitochondria(genome)
        if held:
            raise MitochondrionError(f"{genome} already holds a mitochondrion ({', '.join(held)}); "
                                     f"remove mitochondrion: {accession} from its config entry")
        assembly_check = "the assembly holds no mitochondrion (NCBI sequence report)"
    else:
        assembly_check = "not checked: a local FASTA cannot be asked; make sure it holds no mitochondrion"

    text = fasta(accession)
    lines = text.splitlines()
    headers = [l for l in lines if l.startswith(">")]
    if len(headers) != 1 or headers[0][1:].split()[0] != accession:
        raise MitochondrionError(f"{accession}: NCBI returned {len(headers)} records instead of this one")
    seq = "".join(l.strip() for l in lines if l and not l.startswith(">"))
    if len(seq) != int(s["slen"]):
        raise MitochondrionError(f"{accession}: {len(seq)} bp returned, NCBI states {s['slen']}")
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_fasta.with_suffix(out_fasta.suffix + ".partial")
    body = f">{accession} {s.get('title', '')}\n" + "".join(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
    tmp.write_text(body)
    tmp.rename(out_fasta)
    return {"retrieved": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "source": "ncbi_nucleotide",
            "accession": accession, "title": s.get("title"), "organism": s.get("organism"),
            "tax_id": mt_taxid, "sourcedb": s.get("sourcedb"), "length": len(seq),
            "md5": hashlib.md5(body.encode()).hexdigest(), "genome": genome, "genome_tax_id": g_taxid,
            "species_tax_id": species_of(g_taxid, ranks), "assembly_check": assembly_check,
            "screened": "no: one curated mitochondrial genome of the host's own species"}

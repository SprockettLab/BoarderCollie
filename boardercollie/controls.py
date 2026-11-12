"""The foreign controls (resources/foreign_controls.tsv): known foreign sequences added to
every genome's queries, fetched once into work_dir/foreign_controls.fna."""

from pathlib import Path


def listed(tsv) -> list[str]:
    """Accessions in the controls table, in order."""
    return [line.split("\t")[0] for line in Path(tsv).read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def fetched(fna) -> list[str]:
    """Accessions in a fetched controls FASTA (headers >FOREIGN_<accession>), in order."""
    return [line[len(">FOREIGN_"):].strip() for line in Path(fna).read_text().splitlines()
            if line.startswith(">FOREIGN_")]


def up_to_date(tsv, fna) -> bool:
    """True when the fetched FASTA exists and holds exactly the listed controls.

    The table lives in BoarderCollie's own folder, and every calling workflow (each
    MAGmaker run folder) has its own fresh clone, so the table's mtime is the clone's,
    not a change. The fetched FASTA is shared; when it already holds the listed
    accessions, the table's mtime must not make every downstream step look stale."""
    fna = Path(fna)
    return fna.exists() and fetched(fna) == listed(tsv)

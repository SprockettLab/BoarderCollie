"""Reading and checking the BoarderCollie config file.

Every threshold the workflow uses is a setting with a default, listed in
DEFAULT_SETTINGS. A config may override any of them; an unknown setting name
is an error, so a misspelled setting can never be silently ignored.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .references import is_automatic

# An NCBI Nucleotide accession with its version: NC_012920.1, CM045128.1, PZ236456.1.
NUCLEOTIDE_ACCESSION = re.compile(r"[A-Z]{1,4}_?[0-9]{5,10}\.[0-9]+")

DEFAULT_SETTINGS: dict[str, Any] = {
    # Download the genomes and join them as they are: no screening, no checks, no cleaning.
    "download_only": False,
    # Each check can be turned off. A region a check would have kept goes on to the
    # remaining checks and is removed unless one of them keeps it. fcs_gx off is the
    # same as download_only.
    "checks": {
        "fcs_gx": True,           # screen with FCS-GX (off: nothing flagged, nothing cleaned)
        "retrovirus": True,       # rule 2
        "related_genome": True,   # rule 3
        "lineage_repeat": True,   # rule 4
        "repeat_family": True,    # rule 5
    },
    # Rule 2: FCS-GX hits whose best-match tax-id falls in these NCBI taxa are
    # kept untested. 11632 = Retroviridae (endogenous retroviruses).
    "keep_taxa": [11632],
    # Aligned fraction (segment bases in at least one alignment / length).
    "host_like_min": 0.5,
    "foreign_like_max": 0.1,
    # Rule 3: a match to the closest relative counts only when the relative's
    # sequence extends at least this far beyond the match on one side.
    "embedded_flank_bp": 10000,
    # Rule 4: host-lineage tandem repeat.
    "repeat_trf_min": 0.8,
    "repeat_copy_min_cover": 0.5,
    "repeat_lineage_rank": "genus",
    # Rule 5: repeat family.
    "family_repeat_min": 0.5,
    "family_unit_min_identity": 0.8,
    "family_period_tolerance": 0.1,
    # Regions shorter than this get FCS-GX's own action: too short to align reliably.
    "min_tested_length": 300,
    # Controls.
    "controls_per_segment": 3,
    "background_per_segment": 100,
    "control_pass_min": 0.8,
    "control_max_n_frac": 0.1,
    "seed": 20261104,
    # Rule 1: what happens to a genome whose host controls fail.
    "on_controls_failed": "fcs_only",
    # Automatic choice of comparison genomes (when the config leaves them blank).
    # The closest relative and every comparison genome need at least this scaffold N50.
    "relative_min_scaffold_n50": 1000000,
    # One comparison genome per rank, nearest first; the last is rule 4's distant reference.
    "ladder_ranks": ["family", "order"],
    # Supporting evidence (reported, never deciding).
    "conflict_coding_min": 0.7,
    "conflict_background_percentiles": [1, 99],
}

FRACTION_SETTINGS = (
    "host_like_min", "foreign_like_max", "repeat_trf_min", "repeat_copy_min_cover",
    "family_repeat_min", "family_unit_min_identity", "family_period_tolerance",
    "control_pass_min", "control_max_n_frac", "conflict_coding_min",
)
LINEAGE_RANKS = ("species", "genus", "family", "order")
LADDER_RANK_ORDER = ("genus", "family", "order", "class")
INDEX_TYPES = ("bowtie2", "minimap2")
ON_CONTROLS_FAILED = ("fcs_only", "stop")


class ConfigError(ValueError):
    """The config file is missing something or has an invalid value."""


@dataclass
class GenomeSpec:
    """One genome of a host: an NCBI accession, or a local FASTA with its tax-id."""

    accession: str | None = None
    fasta: str | None = None
    tax_id: int | None = None
    relative: str | None = None
    # The species' complete mitochondrial genome (an NCBI Nucleotide accession, e.g. the
    # RefSeq NC_ record), added to the host reference after this genome. For assemblies
    # deposited without their mitochondrion, as most are: NCBI adds the same record to
    # its RefSeq copies. Checked when fetched (fetch_mitochondrion): it must be a complete
    # mitochondrion of the same species, and the assembly must not already hold one.
    mitochondrion: str | None = None

    @property
    def key(self) -> str:
        return self.accession or Path(self.fasta).stem


@dataclass
class HostSpec:
    name: str
    genomes: list[GenomeSpec]
    ladder: list[str]
    indexes: list[str] = field(default_factory=lambda: list(INDEX_TYPES))
    # A fixed path the joined host reference is also written to (hard link, or copy
    # across file systems), for a workflow that needs the path before the run, such
    # as MAGmaker's host_filter genome.
    reference_out: str | None = None


@dataclass
class Config:
    library_dir: Path
    index_dir: Path
    work_dir: Path
    fcs_gx_db: Path | None
    hosts: list[HostSpec]
    settings: dict[str, Any]
    fcs_gx_db_build: str | None = None
    # A folder of FCS-GX reports made earlier (boardercollie.precomputed); genomes
    # with reports there are not screened again.
    fcs_gx_reports: Path | None = None

    def settings_hash(self) -> str:
        """Short hash of the settings, part of every result's run key."""
        return settings_hash(self.settings)


def settings_hash(settings: dict[str, Any]) -> str:
    text = json.dumps(settings, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def write_if_changed(path: str | Path, text: str) -> bool:
    """Write text to path unless the file already holds exactly that text.

    The text is written to a temporary file in the same folder and renamed over the
    path, so a reader sees either the old file or the new one, never part of one.
    Returns True when the file was written.
    """
    path = Path(path)
    try:
        if path.read_text() == text:
            return False
    except FileNotFoundError:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        # mkstemp makes the file readable by its owner only; give it the usual permissions.
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise
    return True


def merge_settings(overrides: dict[str, Any] | None) -> dict[str, Any]:
    """Defaults with overrides applied; checks names, types and ranges."""
    settings = copy.deepcopy(DEFAULT_SETTINGS)
    for key, value in (overrides or {}).items():
        if key not in DEFAULT_SETTINGS:
            raise ConfigError(f"unknown setting: {key}")
        if key == "checks":
            if not isinstance(value, dict):
                raise ConfigError("checks must map check names to true or false")
            unknown = set(value) - set(DEFAULT_SETTINGS["checks"])
            if unknown:
                raise ConfigError(f"unknown checks: {sorted(unknown)}; known: {sorted(DEFAULT_SETTINGS['checks'])}")
            settings["checks"].update(value)
        else:
            settings[key] = value
    if not isinstance(settings["download_only"], bool) or not all(isinstance(v, bool) for v in settings["checks"].values()):
        raise ConfigError("download_only and every check must be true or false")
    if settings["download_only"] or not settings["checks"]["fcs_gx"]:
        settings["download_only"] = True
        settings["checks"] = {k: False for k in settings["checks"]}
    for key in FRACTION_SETTINGS:
        v = settings[key]
        if not isinstance(v, (int, float)) or not 0 <= v <= 1:
            raise ConfigError(f"{key} must be a number from 0 to 1, got {v!r}")
    if settings["foreign_like_max"] >= settings["host_like_min"]:
        raise ConfigError("foreign_like_max must be below host_like_min")
    if settings["repeat_lineage_rank"] not in LINEAGE_RANKS:
        raise ConfigError(f"repeat_lineage_rank must be one of {LINEAGE_RANKS}")
    if settings["on_controls_failed"] not in ON_CONTROLS_FAILED:
        raise ConfigError(f"on_controls_failed must be one of {ON_CONTROLS_FAILED}")
    for key in ("embedded_flank_bp", "controls_per_segment", "background_per_segment", "seed", "min_tested_length"):
        if not isinstance(settings[key], int) or settings[key] < 0:
            raise ConfigError(f"{key} must be a non-negative integer")
    if not settings["controls_per_segment"] or not settings["background_per_segment"]:
        raise ConfigError("controls_per_segment and background_per_segment must be at least 1")
    if not all(isinstance(t, int) for t in settings["keep_taxa"]):
        raise ConfigError("keep_taxa must be a list of NCBI tax-ids (integers)")
    if not isinstance(settings["relative_min_scaffold_n50"], int) or settings["relative_min_scaffold_n50"] < 0:
        raise ConfigError("relative_min_scaffold_n50 must be a non-negative integer")
    lr = settings["ladder_ranks"]
    if (not isinstance(lr, list) or not lr or not set(lr) <= set(LADDER_RANK_ORDER) or len(set(lr)) != len(lr)
            or lr != sorted(lr, key=LADDER_RANK_ORDER.index)):
        raise ConfigError(f"ladder_ranks must list ranks from {LADDER_RANK_ORDER}, nearest first, each once")
    lo, hi = settings["conflict_background_percentiles"]
    if not 0 <= lo < hi <= 100:
        raise ConfigError("conflict_background_percentiles must be [low, high] with 0 <= low < high <= 100")
    return settings


def _genome(entry: Any, host: str) -> GenomeSpec:
    if isinstance(entry, str):
        return GenomeSpec(accession=entry)
    if not isinstance(entry, dict):
        raise ConfigError(f"host {host}: each genome is an accession or a mapping")
    unknown = set(entry) - {"accession", "fasta", "tax_id", "relative", "mitochondrion"}
    if unknown:
        raise ConfigError(f"host {host}: unknown genome keys {sorted(unknown)}")
    g = GenomeSpec(**entry)
    if is_automatic(g.relative):
        g.relative = None
    if bool(g.accession) == bool(g.fasta):
        raise ConfigError(f"host {host}: give either accession or fasta for each genome")
    if g.fasta and not isinstance(g.tax_id, int):
        raise ConfigError(f"host {host}: a local fasta needs an integer tax_id")
    if g.mitochondrion is not None and not (isinstance(g.mitochondrion, str)
                                            and NUCLEOTIDE_ACCESSION.fullmatch(g.mitochondrion)):
        raise ConfigError(f"host {host}: genome {g.key}: mitochondrion must be an NCBI Nucleotide accession "
                          f"with its version, e.g. NC_012920.1; got {g.mitochondrion!r}")
    return g


def parse_config(data: dict[str, Any]) -> Config:
    settings = merge_settings(data.get("settings"))
    for key in ("library_dir", "index_dir", "work_dir", "hosts"):
        if key not in data:
            raise ConfigError(f"missing required key: {key}")
    # The FCS-GX database is needed unless nothing is screened or reports made earlier
    # are given; whether those reports cover every genome is checked when the workflow starts.
    if not settings["download_only"] and not data.get("fcs_gx_db") and not data.get("fcs_gx_reports"):
        raise ConfigError("missing required key: fcs_gx_db (or fcs_gx_reports, a folder of FCS-GX reports made earlier)")
    unknown = set(data) - {"library_dir", "index_dir", "work_dir", "fcs_gx_db", "fcs_gx_db_build", "fcs_gx_reports",
                           "hosts", "settings"}
    if unknown:
        raise ConfigError(f"unknown top-level keys: {sorted(unknown)}")
    if not isinstance(data["hosts"], list) or not data["hosts"]:
        raise ConfigError("hosts must be a non-empty list")
    hosts, names, outs = [], set(), set()
    for h in data["hosts"]:
        name = h.get("name") if isinstance(h, dict) else None
        if not name:
            raise ConfigError("every host needs a name")
        if name in names:
            raise ConfigError(f"duplicate host name: {name}")
        names.add(name)
        genomes = [_genome(g, name) for g in h.get("genomes") or []]
        if not genomes:
            raise ConfigError(f"host {name}: genomes must list at least one genome")
        keys = [g.key for g in genomes]
        if len(set(keys)) != len(keys):
            raise ConfigError(f"host {name}: a genome is listed twice")
        indexes = h.get("indexes", list(INDEX_TYPES))
        if not set(indexes) <= set(INDEX_TYPES):
            raise ConfigError(f"host {name}: indexes must be from {INDEX_TYPES}")
        # A blank, none, NA or auto relative or ladder is chosen automatically
        # from the NCBI taxonomy when the workflow starts (boardercollie.references).
        ladder = h.get("ladder")
        if is_automatic(ladder) or ladder == [] or (isinstance(ladder, list) and all(is_automatic(a) for a in ladder)):
            ladder = []
        elif not isinstance(ladder, list) or not all(isinstance(a, str) and not is_automatic(a) for a in ladder):
            raise ConfigError(f"host {name}: ladder must list reference accessions, nearest first, or be left blank")
        for g in genomes:
            if g.relative is None:
                continue
            if g.relative == g.key or g.relative in ladder:
                raise ConfigError(f"host {name}: genome {g.key}: the relative must differ from the genome and the ladder")
        if set(ladder) & set(keys):
            raise ConfigError(f"host {name}: a ladder reference is also one of the host's genomes")
        unknown = set(h) - {"name", "genomes", "ladder", "indexes", "reference_out"}
        if unknown:
            raise ConfigError(f"host {name}: unknown keys {sorted(unknown)}")
        out = h.get("reference_out")
        if out is not None and (not isinstance(out, str) or not out.strip()):
            raise ConfigError(f"host {name}: reference_out must be a file path")
        if out and out in outs:
            raise ConfigError(f"host {name}: reference_out {out} is also another host's")
        outs.add(out)
        hosts.append(HostSpec(name=name, genomes=genomes, ladder=list(ladder), indexes=list(indexes),
                              reference_out=out or None))
    return Config(
        library_dir=Path(data["library_dir"]),
        index_dir=Path(data["index_dir"]),
        work_dir=Path(data["work_dir"]),
        fcs_gx_db=Path(data["fcs_gx_db"]) if data.get("fcs_gx_db") else None,
        hosts=hosts,
        settings=settings,
        fcs_gx_db_build=str(data["fcs_gx_db_build"]) if data.get("fcs_gx_db_build") else None,
        fcs_gx_reports=Path(data["fcs_gx_reports"]) if data.get("fcs_gx_reports") else None,
    )


def load_config(path: str | Path) -> Config:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: not a YAML mapping")
    return parse_config(data)

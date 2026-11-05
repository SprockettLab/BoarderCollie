"""Reading and checking the BoarderCollie config file.

Every threshold the workflow uses is a setting with a default, listed in
DEFAULT_SETTINGS. A config may override any of them; an unknown setting name
is an error, so a misspelled setting can never be silently ignored.
"""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_SETTINGS: dict[str, Any] = {
    # Rule R0: FCS-GX hits whose best-match tax-id falls in these NCBI taxa are
    # kept untested. 11632 = Retroviridae (endogenous retroviruses).
    "keep_taxa": [11632],
    # Aligned fraction (segment bases in at least one alignment / length).
    "host_like_min": 0.5,
    "foreign_like_max": 0.1,
    # Rule R1: a match to the closest relative counts only when the relative's
    # sequence extends at least this far beyond the match on one side.
    "embedded_flank_bp": 10000,
    # Rule R2: host-lineage tandem repeat.
    "repeat_trf_min": 0.8,
    "repeat_copy_min_cover": 0.5,
    "repeat_lineage_rank": "genus",
    # Rule R2b: repeat family.
    "family_repeat_min": 0.5,
    "family_unit_min_identity": 0.8,
    "family_period_tolerance": 0.1,
    # Controls.
    "controls_per_segment": 3,
    "background_per_segment": 100,
    "control_pass_min": 0.8,
    "control_max_n_frac": 0.1,
    "seed": 20261104,
    # Rule C: what happens to a genome whose host controls fail.
    "on_controls_failed": "fcs_only",
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

    @property
    def key(self) -> str:
        return self.accession or Path(self.fasta).stem


@dataclass
class HostSpec:
    name: str
    genomes: list[GenomeSpec]
    ladder: list[str]
    indexes: list[str] = field(default_factory=lambda: list(INDEX_TYPES))


@dataclass
class Config:
    library_dir: Path
    index_dir: Path
    work_dir: Path
    fcs_gx_db: Path
    hosts: list[HostSpec]
    settings: dict[str, Any]
    fcs_gx_stage: Path | None = None
    fcs_gx_db_build: str | None = None

    def settings_hash(self) -> str:
        """Short hash of the settings, part of every result's run key."""
        return settings_hash(self.settings)


def settings_hash(settings: dict[str, Any]) -> str:
    text = json.dumps(settings, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def merge_settings(overrides: dict[str, Any] | None) -> dict[str, Any]:
    """Defaults with overrides applied; checks names, types and ranges."""
    settings = copy.deepcopy(DEFAULT_SETTINGS)
    for key, value in (overrides or {}).items():
        if key not in DEFAULT_SETTINGS:
            raise ConfigError(f"unknown setting: {key}")
        settings[key] = value
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
    for key in ("embedded_flank_bp", "controls_per_segment", "background_per_segment", "seed"):
        if not isinstance(settings[key], int) or settings[key] < 0:
            raise ConfigError(f"{key} must be a non-negative integer")
    if not settings["controls_per_segment"] or not settings["background_per_segment"]:
        raise ConfigError("controls_per_segment and background_per_segment must be at least 1")
    if not all(isinstance(t, int) for t in settings["keep_taxa"]):
        raise ConfigError("keep_taxa must be a list of NCBI tax-ids (integers)")
    lo, hi = settings["conflict_background_percentiles"]
    if not 0 <= lo < hi <= 100:
        raise ConfigError("conflict_background_percentiles must be [low, high] with 0 <= low < high <= 100")
    return settings


def _genome(entry: Any, host: str) -> GenomeSpec:
    if isinstance(entry, str):
        return GenomeSpec(accession=entry)
    if not isinstance(entry, dict):
        raise ConfigError(f"host {host}: each genome is an accession or a mapping")
    unknown = set(entry) - {"accession", "fasta", "tax_id", "relative"}
    if unknown:
        raise ConfigError(f"host {host}: unknown genome keys {sorted(unknown)}")
    g = GenomeSpec(**entry)
    if bool(g.accession) == bool(g.fasta):
        raise ConfigError(f"host {host}: give either accession or fasta for each genome")
    if g.fasta and not isinstance(g.tax_id, int):
        raise ConfigError(f"host {host}: a local fasta needs an integer tax_id")
    return g


def parse_config(data: dict[str, Any]) -> Config:
    for key in ("library_dir", "index_dir", "work_dir", "fcs_gx_db", "hosts"):
        if key not in data:
            raise ConfigError(f"missing required key: {key}")
    unknown = set(data) - {"library_dir", "index_dir", "work_dir", "fcs_gx_db", "fcs_gx_stage", "fcs_gx_db_build", "hosts", "settings"}
    if unknown:
        raise ConfigError(f"unknown top-level keys: {sorted(unknown)}")
    if not isinstance(data["hosts"], list) or not data["hosts"]:
        raise ConfigError("hosts must be a non-empty list")
    hosts, names = [], set()
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
        # Automatic choice of relatives and ladder references from the NCBI
        # taxonomy is not implemented yet; until it is, both are given.
        ladder = h.get("ladder") or []
        if not ladder or not all(isinstance(a, str) for a in ladder):
            raise ConfigError(f"host {name}: ladder must list reference accessions, nearest first")
        for g in genomes:
            if not g.relative:
                raise ConfigError(f"host {name}: genome {g.key} needs a relative (automatic choice not available yet)")
            if g.relative == g.key or g.relative in ladder:
                raise ConfigError(f"host {name}: genome {g.key}: the relative must differ from the genome and the ladder")
        if set(ladder) & set(keys):
            raise ConfigError(f"host {name}: a ladder reference is also one of the host's genomes")
        hosts.append(HostSpec(name=name, genomes=genomes, ladder=list(ladder), indexes=list(indexes)))
    return Config(
        library_dir=Path(data["library_dir"]),
        index_dir=Path(data["index_dir"]),
        work_dir=Path(data["work_dir"]),
        fcs_gx_db=Path(data["fcs_gx_db"]),
        hosts=hosts,
        settings=merge_settings(data.get("settings")),
        fcs_gx_stage=Path(data["fcs_gx_stage"]) if data.get("fcs_gx_stage") else None,
        fcs_gx_db_build=str(data["fcs_gx_db_build"]) if data.get("fcs_gx_db_build") else None,
    )


def load_config(path: str | Path) -> Config:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: not a YAML mapping")
    return parse_config(data)

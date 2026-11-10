"""The decision rules: which flagged regions of a host genome are removed.

After FCS-GX, every flagged region goes through the relative-genome alignment
test: does it align to related host genomes (the closest relative and the
reference ladder) as well as the same genome's own host DNA does? The rules,
applied in this order, are:

  1  Control check, per genome. At least control_pass_min of the host
     controls must be host_like and every foreign control foreign_like.
     If not, the relative-genome alignment test cannot be applied to this
     genome; on_controls_failed decides what happens (fcs_only or stop).
  2  Endogenous retrovirus: FCS-GX's best-match tax-id lies within a taxon in
     keep_taxa. Kept, untested.
  3  Host by alignment: aligned fraction >= host_like_min to a ladder
     reference, or to the closest relative with the match embedded in a
     longer relative sequence. Kept.
  4  Host-lineage tandem repeat: TRF covers >= repeat_trf_min of the region,
     at least one copy in another assembly of the same lineage, no copy in
     the most distant reference. Kept.
  5  Repeat family: the region's main tandem repeat covers
     >= family_repeat_min of it and its unit matches a unit confirmed by
     rule 3 or 4 in the same genome. Kept.
  6  Everything else. Removed.

Regions shorter than min_tested_length are too short for the alignment test
(whole-genome alignment does not align short sequences reliably: in the first
test runs, under half of 200-300 bp host controls aligned, against 99-100% of
longer ones). After rule 2, such a region gets FCS-GX's own action.

Only the host controls of regions that go through the test count in rule 1;
controls drawn for regions kept by rule 2 or too short do not.

The rule column of decisions.tsv holds the number of the rule that decided
each region, "too_short" for a region under min_tested_length, or "fcs_only"
when rule 1 failed and FCS-GX's own action was applied.

Composition and gene evidence never decide. It is compared with the decision
and any disagreement is reported in evidence_conflict.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable

from .repeats import unit_identity

TEST_APPLIED = "applied"
TEST_NOT_APPLIED = "not applied: host controls failed"
TEST_NOT_NEEDED = "not needed: FCS-GX flagged no regions"
TEST_NONE_TESTABLE = "not needed: every flagged region kept as a retrovirus or too short for the test"
TEST_NO_RELATIVE = "not applied: no relative genome found"
TEST_NO_LADDER = "not applied: no comparison genomes found"
FCS_REMOVING_ACTIONS = ("EXCLUDE", "TRIM", "FIX")


@dataclass
class Alignment:
    """One region (or control) aligned to one reference genome."""

    reference: str          # accession of the reference genome
    role: str               # "relative" or "ladder"
    aligned_fraction: float
    identity: float = 0.0
    embedded: bool = False  # relative match lies inside a longer relative sequence


@dataclass
class RepeatEvidence:
    trf_cover: float           # fraction of the region in any tandem repeat
    main_repeat_cover: float   # fraction covered by the repeat covering the most bases
    main_unit: str             # consensus unit of that repeat ("" if none)
    lineage_copies: int        # copies in other assemblies of the same lineage
    distant_copies: int        # copies in the most distant ladder reference


@dataclass
class CompositionEvidence:
    cpg_oe: float
    cpg_oe_background_low: float    # background percentiles (conflict_background_percentiles)
    cpg_oe_background_high: float
    coding_fraction: float | None = None     # Prodigal; None until measured
    prokaryotic_rrna: bool | None = None     # barrnap; None until measured


@dataclass
class Region:
    region_id: str
    accession: str
    seq_id: str
    start: int
    end: int
    fcs_action: str
    fcs_div: str
    fcs_top_tax: str
    fcs_top_taxid: int | None
    source: str
    alignments: list[Alignment] = field(default_factory=list)
    repeats: RepeatEvidence | None = None
    composition: CompositionEvidence | None = None

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Control:
    control_id: str
    kind: str                # "host" (from this genome) or "foreign"
    alignments: list[Alignment] = field(default_factory=list)
    region_id: str | None = None  # the flagged region a host control was drawn for


@dataclass
class Decision:
    region: Region
    best_aligned_fraction: float
    align_class: str
    rule: str
    conclusion: str           # "host" or "foreign"
    action: str               # "keep" or "remove"
    relative_alignment_test: str
    evidence_conflict: bool = False
    conflict_reasons: str = ""
    note: str = ""


@dataclass
class GenomeResult:
    accession: str
    status: str               # "decided", "fcs_only" or "stopped"
    relative_alignment_test: str
    host_controls_passed: int
    host_controls_total: int
    foreign_controls_passed: int
    foreign_controls_total: int
    decisions: list[Decision]


def best_fraction(alignments: Iterable[Alignment]) -> float:
    return max((a.aligned_fraction for a in alignments), default=0.0)


def align_class(fraction: float, settings: dict) -> str:
    if fraction >= settings["host_like_min"]:
        return "host_like"
    if fraction < settings["foreign_like_max"]:
        return "foreign_like"
    return "partial"


def control_check(controls: list[Control], settings: dict) -> tuple[bool, int, int, int, int]:
    host = [c for c in controls if c.kind == "host"]
    foreign = [c for c in controls if c.kind == "foreign"]
    host_ok = sum(align_class(best_fraction(c.alignments), settings) == "host_like" for c in host)
    foreign_ok = sum(align_class(best_fraction(c.alignments), settings) == "foreign_like" for c in foreign)
    passed = (bool(host) and host_ok / len(host) >= settings["control_pass_min"]
              and foreign_ok == len(foreign))
    return passed, host_ok, len(host), foreign_ok, len(foreign)


def is_kept_taxon(taxid: int | None, settings: dict, lineage: Callable[[int], set[int]]) -> bool:
    """True if taxid is, or descends from, a taxon in keep_taxa (NCBI taxonomy)."""
    if taxid is None:
        return False
    return bool(({taxid} | lineage(taxid)) & set(settings["keep_taxa"]))


def _host_by_alignment(r: Region, settings: dict) -> bool:
    for a in r.alignments:
        if a.aligned_fraction < settings["host_like_min"]:
            continue
        if a.role == "ladder" or (a.role == "relative" and a.embedded):
            return True
    return False


def _host_lineage_repeat(r: Region, settings: dict) -> bool:
    x = r.repeats
    return (x is not None and x.trf_cover >= settings["repeat_trf_min"]
            and x.lineage_copies >= 1 and x.distant_copies == 0)


def _conflict(r: Region, action: str, settings: dict) -> tuple[bool, str]:
    c = r.composition
    if c is None:
        return False, ""
    foreign_signs = []
    if c.cpg_oe > c.cpg_oe_background_high:
        foreign_signs.append(f"CpG o/e {c.cpg_oe:.3f} above host background ({c.cpg_oe_background_high:.3f})")
    if c.coding_fraction is not None and c.coding_fraction >= settings["conflict_coding_min"]:
        foreign_signs.append(f"coding fraction {c.coding_fraction:.2f}")
    if c.prokaryotic_rrna:
        foreign_signs.append("prokaryotic rRNA gene")
    if action == "keep" and foreign_signs:
        return True, "kept, but " + "; ".join(foreign_signs)
    if action == "remove" and not foreign_signs and c.cpg_oe >= c.cpg_oe_background_low:
        return True, "removed, but composition and genes are within the host background"
    return False, ""


def _decision(r: Region, rule: str, conclusion: str, action: str, test: str, settings: dict, note: str = "") -> Decision:
    frac = best_fraction(r.alignments)
    conflict, reasons = _conflict(r, action, settings)
    return Decision(r, frac, align_class(frac, settings), rule, conclusion, action, test, conflict, reasons, note)


def _fcs_decision(r: Region, rule: str, test: str, settings: dict, why: str) -> Decision:
    """FCS-GX's own action for a region the alignment test does not decide."""
    if r.fcs_action in FCS_REMOVING_ACTIONS:
        return _decision(r, rule, "foreign", "remove", test, settings, note=f"{why}: FCS-GX action {r.fcs_action} applied")
    return _decision(r, rule, "host", "keep", test, settings, note=f"{why}: FCS-GX action {r.fcs_action}, no removal")


def decide_genome(
    accession: str,
    regions: list[Region],
    controls: list[Control],
    settings: dict,
    lineage: Callable[[int], set[int]],
    test_unavailable: str | None = None,
) -> GenomeResult:
    """Apply rules 1 to 6 to every flagged region of one genome.

    test_unavailable: why the relative-genome alignment test cannot be run at all
    (TEST_NO_RELATIVE or TEST_NO_LADDER); the genome then gets FCS-GX's own actions,
    as when its host controls fail.
    """
    kept_taxon = {r.region_id for r in regions if is_kept_taxon(r.fcs_top_taxid, settings, lineage)}
    too_short = {r.region_id for r in regions
                 if r.region_id not in kept_taxon and r.length < settings["min_tested_length"]}
    tested = {r.region_id for r in regions} - kept_taxon - too_short
    counted = [c for c in controls if c.kind != "host" or c.region_id is None or c.region_id in tested]
    passed, h_ok, h_n, f_ok, f_n = control_check(counted, settings)
    if not regions:
        # Nothing flagged: there is nothing to test and no host controls are drawn.
        return GenomeResult(accession, "no flagged regions", TEST_NOT_NEEDED, h_ok, h_n, f_ok, f_n, [])
    short_label = f"not applied: shorter than {settings['min_tested_length']} bp"
    if not tested:
        out = []
        for r in regions:
            if r.region_id in kept_taxon:
                out.append(_decision(r, "2", "host", "keep", TEST_NONE_TESTABLE, settings))
            else:
                out.append(_fcs_decision(r, "too_short", short_label, settings, "too short for the alignment test"))
        return GenomeResult(accession, "decided", TEST_NONE_TESTABLE, h_ok, h_n, f_ok, f_n, out)

    if test_unavailable or not passed:
        label = test_unavailable or TEST_NOT_APPLIED
        if settings["on_controls_failed"] == "stop":
            return GenomeResult(accession, "stopped", label, h_ok, h_n, f_ok, f_n, [])
        out = []
        for r in regions:
            if is_kept_taxon(r.fcs_top_taxid, settings, lineage):
                out.append(_decision(r, "2", "host", "keep", label, settings))
            elif r.fcs_action in FCS_REMOVING_ACTIONS:
                out.append(_decision(r, "fcs_only", "foreign", "remove", label, settings,
                                     note=f"FCS-GX action {r.fcs_action} applied without the relative-genome alignment test"))
            else:
                out.append(_decision(r, "fcs_only", "host", "keep", label, settings,
                                     note=f"FCS-GX action {r.fcs_action}: no removal"))
        return GenomeResult(accession, "fcs_only", label, h_ok, h_n, f_ok, f_n, out)

    decided: dict[str, Decision] = {}
    pending: list[Region] = []
    for r in regions:
        if r.region_id in kept_taxon:
            decided[r.region_id] = _decision(r, "2", "host", "keep", TEST_APPLIED, settings)
        elif r.region_id in too_short:
            decided[r.region_id] = _fcs_decision(r, "too_short", short_label, settings, "too short for the alignment test")
        elif _host_by_alignment(r, settings):
            decided[r.region_id] = _decision(r, "3", "host", "keep", TEST_APPLIED, settings)
        elif _host_lineage_repeat(r, settings):
            decided[r.region_id] = _decision(r, "4", "host", "keep", TEST_APPLIED, settings)
        else:
            pending.append(r)

    # Rule 5: only units confirmed by rule 3 or 4 count, so family members cannot vouch for each other.
    confirmed = [d.region.repeats.main_unit for d in decided.values()
                 if d.rule in ("3", "4") and d.region.repeats and d.region.repeats.main_unit]
    for r in pending:
        x = r.repeats
        match = None
        if x and x.main_unit and x.main_repeat_cover >= settings["family_repeat_min"] and confirmed:
            best = max(((u, unit_identity(x.main_unit, u, settings["family_period_tolerance"])) for u in confirmed),
                       key=lambda t: t[1])
            if best[1] >= settings["family_unit_min_identity"]:
                match = best
        if match:
            decided[r.region_id] = _decision(r, "5", "host", "keep", TEST_APPLIED, settings,
                                             note=f"unit identity {match[1]:.3f} to confirmed unit of {len(match[0])} bp")
        else:
            decided[r.region_id] = _decision(r, "6", "foreign", "remove", TEST_APPLIED, settings)

    return GenomeResult(accession, "decided", TEST_APPLIED, h_ok, h_n, f_ok, f_n,
                        [decided[r.region_id] for r in regions])

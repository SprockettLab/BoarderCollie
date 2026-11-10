"""Decision rules on small made-up evidence, one test per rule and per edge case."""

import pytest

from boardercollie.config import merge_settings
from boardercollie.rules import (
    TEST_APPLIED, TEST_NOT_APPLIED, TEST_NOT_NEEDED, Alignment, CompositionEvidence, Control, Region, RepeatEvidence, decide_genome,
)

RETROVIRIDAE = 11632
LINEAGES = {39068: {RETROVIRIDAE, 10239, 1}, 816: {976, 2, 131567, 1}}


def lineage(taxid):
    return LINEAGES.get(taxid, set())


def region(rid, alignments=(), repeats=None, composition=None, action="EXCLUDE", taxid=816, tax="Bacteroides"):
    return Region(rid, "GCA_TEST.1", f"seq_{rid}", 1, 1000, action, "prok:CFB group bacteria", tax, taxid,
                  "action_report", list(alignments), repeats, composition)


def ladder(frac):
    return Alignment("GCA_LADDER.1", "ladder", frac, 0.9)


def relative(frac, embedded):
    return Alignment("GCA_REL.1", "relative", frac, 0.95, embedded)


def good_controls(n=10):
    host = [Control(f"c{i}", "host", [ladder(0.95)]) for i in range(n)]
    return host + [Control("bact", "foreign", [ladder(0.0)])]


@pytest.fixture
def settings():
    return merge_settings({})


def decide(regions, settings, controls=None):
    return decide_genome("GCA_TEST.1", regions, controls or good_controls(), settings, lineage)


def by_id(result):
    return {d.region.region_id: d for d in result.decisions}


def test_foreign_region_is_removed(settings):
    d = by_id(decide([region("a", [ladder(0.0), relative(0.0, False)])], settings))["a"]
    assert (d.rule, d.action, d.align_class) == ("6", "remove", "foreign_like")
    assert d.relative_alignment_test == TEST_APPLIED


def test_partial_region_is_removed(settings):
    d = by_id(decide([region("a", [ladder(0.3)])], settings))["a"]
    assert (d.rule, d.action, d.align_class) == ("6", "remove", "partial")


def test_ladder_alignment_keeps(settings):
    d = by_id(decide([region("a", [ladder(0.75)])], settings))["a"]
    assert (d.rule, d.action) == ("3", "keep")


def test_relative_alignment_keeps_only_when_embedded(settings):
    res = by_id(decide([region("emb", [relative(1.0, True)]), region("short", [relative(1.0, False)])], settings))
    assert (res["emb"].rule, res["emb"].action) == ("3", "keep")
    assert (res["short"].rule, res["short"].action) == ("6", "remove")


def test_retrovirus_kept_by_taxonomy_not_name(settings):
    d = by_id(decide([region("erv", [ladder(0.0)], taxid=39068, tax="Simian retrovirus 2")], settings))["erv"]
    assert (d.rule, d.action) == ("2", "keep")
    # A name that says "retrovirus" but a tax-id outside Retroviridae is not kept by rule 2.
    d = by_id(decide([region("fake", [ladder(0.0)], taxid=816, tax="not a retrovirus")], settings))["fake"]
    assert d.rule == "6"


def test_host_lineage_repeat_keeps(settings):
    rep = RepeatEvidence(0.95, 0.95, "ACGTTGCAACGTTGCA", lineage_copies=2, distant_copies=0)
    d = by_id(decide([region("rep", [relative(1.0, False)], repeats=rep)], settings))["rep"]
    assert (d.rule, d.action) == ("4", "keep")


def test_repeat_with_copy_in_distant_reference_is_not_r2(settings):
    rep = RepeatEvidence(0.95, 0.95, "ACGT", lineage_copies=2, distant_copies=1)
    d = by_id(decide([region("rep", [ladder(0.0)], repeats=rep)], settings))["rep"]
    assert d.rule == "6"


def test_repeat_family_rotation_and_reverse_complement(settings):
    unit = "ACGTTGCAAGGCTTAC"
    rotated_rc = "GTAAGCCTTGCAACGT"[5:] + "GTAAGCCTTGCAACGT"[:5]
    confirmed = RepeatEvidence(1.0, 1.0, unit, lineage_copies=1, distant_copies=0)
    member = RepeatEvidence(1.0, 1.0, rotated_rc, lineage_copies=0, distant_copies=0)
    res = by_id(decide([region("conf", [ladder(0.0)], repeats=confirmed),
                        region("mem", [ladder(0.0)], repeats=member)], settings))
    assert res["conf"].rule == "4"
    assert (res["mem"].rule, res["mem"].action) == ("5", "keep")


def test_repeat_family_needs_main_repeat_cover(settings):
    confirmed = RepeatEvidence(1.0, 1.0, "ACGTTGCAAGGCTTAC", 1, 0)
    weak = RepeatEvidence(0.3, 0.3, "ACGTTGCAAGGCTTAC", 0, 0)
    res = by_id(decide([region("conf", repeats=confirmed), region("weak", repeats=weak)], settings))
    assert res["weak"].rule == "6"


def test_repeat_family_members_cannot_vouch_for_each_other(settings):
    # Two identical units, neither confirmed by rule 3 or 4: both removed.
    a = RepeatEvidence(1.0, 1.0, "ACGTTGCAAGGCTTAC", 0, 0)
    res = by_id(decide([region("x", repeats=a), region("y", repeats=a)], settings))
    assert res["x"].rule == res["y"].rule == "6"


def test_controls_failed_fcs_only_is_labelled(settings):
    bad = [Control(f"c{i}", "host", [ladder(0.0)]) for i in range(10)]
    result = decide([region("a", [ladder(0.9)], action="EXCLUDE"),
                     region("b", [ladder(0.0)], action="REVIEW"),
                     region("erv", taxid=39068, action="EXCLUDE")], settings, controls=bad)
    assert result.status == "fcs_only"
    assert result.relative_alignment_test == TEST_NOT_APPLIED
    d = by_id(result)
    assert (d["a"].rule, d["a"].action) == ("fcs_only", "remove")
    assert (d["b"].rule, d["b"].action) == ("fcs_only", "keep")
    assert (d["erv"].rule, d["erv"].action) == ("2", "keep")
    assert all(x.relative_alignment_test == TEST_NOT_APPLIED for x in result.decisions)


def test_controls_failed_stop(settings):
    s = merge_settings({"on_controls_failed": "stop"})
    result = decide([region("a")], s, controls=[Control("c", "host", [ladder(0.0)])])
    assert result.status == "stopped" and result.decisions == []


def test_foreign_control_must_not_align(settings):
    controls = good_controls()
    controls[-1] = Control("bact", "foreign", [ladder(0.6)])
    assert decide([region("a")], settings, controls=controls).status == "fcs_only"


def test_control_pass_threshold(settings):
    # 8 of 10 host controls host_like passes at the default 0.8; 7 of 10 does not.
    def ctl(n_ok):
        return ([Control(f"o{i}", "host", [ladder(0.9)]) for i in range(n_ok)]
                + [Control(f"b{i}", "host", [ladder(0.0)]) for i in range(10 - n_ok)]
                + [Control("bact", "foreign", [])])
    assert decide([region("a")], settings, controls=ctl(8)).status == "decided"
    assert decide([region("a")], settings, controls=ctl(7)).status == "fcs_only"


def test_evidence_conflict_reported_but_not_deciding(settings):
    comp = CompositionEvidence(cpg_oe=1.1, cpg_oe_background_low=0.1, cpg_oe_background_high=0.35,
                               coding_fraction=0.9, prokaryotic_rrna=True)
    d = by_id(decide([region("a", [ladder(0.9)], composition=comp)], settings))["a"]
    assert d.action == "keep" and d.evidence_conflict
    assert "CpG" in d.conflict_reasons and "rRNA" in d.conflict_reasons


def test_thresholds_are_settings(settings):
    s = merge_settings({"host_like_min": 0.8})
    d = by_id(decide([region("a", [ladder(0.75)])], s))["a"]
    assert d.rule == "6"


def test_genome_with_nothing_flagged(settings):
    # No regions, so no host controls: this is not a failed control check.
    result = decide([], settings, controls=[Control("bact", "foreign", [ladder(0.0)])])
    assert result.status == "no flagged regions"
    assert result.relative_alignment_test == TEST_NOT_NEEDED and result.decisions == []


def short_region(rid, length, **kw):
    r = region(rid, **kw)
    r.end = r.start + length - 1
    return r


def test_only_controls_of_tested_regions_count(settings):
    # 30 failing controls belong to a retrovirus region kept untested; the tested region's 3 controls align.
    controls = [Control(f"erv_C{i}", "host", [ladder(0.0)], "erv") for i in range(30)]
    controls += [Control(f"t_C{i}", "host", [ladder(0.95)], "t") for i in range(3)]
    controls += [Control("bact", "foreign", [ladder(0.0)])]
    result = decide([region("t", [ladder(0.9)]), region("erv", taxid=39068)], settings, controls=controls)
    assert result.status == "decided" and result.relative_alignment_test == TEST_APPLIED
    assert (result.host_controls_passed, result.host_controls_total) == (3, 3)
    d = by_id(result)
    assert (d["t"].rule, d["t"].action) == ("3", "keep")
    assert (d["erv"].rule, d["erv"].action) == ("2", "keep")


def test_short_region_gets_fcs_action_and_its_controls_do_not_count(settings):
    controls = [Control(f"s_C{i}", "host", [ladder(0.0)], "s") for i in range(3)]
    controls += [Control(f"t_C{i}", "host", [ladder(0.95)], "t") for i in range(3)]
    controls += [Control("bact", "foreign", [ladder(0.0)])]
    regions = [region("t", [ladder(0.9)]),
               short_region("s", 299, alignments=[ladder(0.0)], action="EXCLUDE"),
               short_region("r", 250, alignments=[ladder(0.0)], action="REVIEW"),
               short_region("ok", 300, alignments=[ladder(0.0)])]
    result = decide(regions, settings, controls=controls)
    assert result.status == "decided" and (result.host_controls_passed, result.host_controls_total) == (3, 3)
    d = by_id(result)
    assert (d["s"].rule, d["s"].action) == ("too_short", "remove")
    assert (d["r"].rule, d["r"].action) == ("too_short", "keep")
    assert d["s"].relative_alignment_test == "not applied: shorter than 300 bp"
    assert "too short" in d["s"].note
    assert (d["ok"].rule, d["ok"].action) == ("6", "remove")  # 300 bp is tested


def test_short_length_is_a_setting(settings):
    s = merge_settings({"min_tested_length": 0})
    result = decide([short_region("s", 100, alignments=[ladder(0.9)])], s)
    assert by_id(result)["s"].rule == "3"


def test_nothing_long_enough_to_test(settings):
    from boardercollie.rules import TEST_NONE_TESTABLE
    result = decide([short_region("s", 200, action="EXCLUDE"), region("erv", taxid=39068)], settings,
                    controls=[Control("bact", "foreign", [ladder(0.0)])])
    assert result.status == "decided" and result.relative_alignment_test == TEST_NONE_TESTABLE
    d = by_id(result)
    assert (d["s"].rule, d["s"].action) == ("too_short", "remove")
    assert (d["erv"].rule, d["erv"].action) == ("2", "keep")

@pytest.mark.parametrize("why", ["not applied: no relative genome found", "not applied: no comparison genomes found"])
def test_no_references_falls_back_to_fcs_only(settings, why):
    from boardercollie.rules import TEST_NO_LADDER, TEST_NO_RELATIVE
    assert why in (TEST_NO_RELATIVE, TEST_NO_LADDER)
    # Controls pass and the region aligns, but without references the test is not run.
    result = decide_genome("GCA_TEST.1", [region("a", [ladder(0.9)], action="EXCLUDE"),
                                          region("erv", taxid=39068, action="EXCLUDE")],
                           good_controls(), settings, lineage, why)
    assert result.status == "fcs_only" and result.relative_alignment_test == why
    d = by_id(result)
    assert (d["a"].rule, d["a"].action) == ("fcs_only", "remove")
    assert (d["erv"].rule, d["erv"].action) == ("2", "keep")
    assert all(x.relative_alignment_test == why for x in result.decisions)


def checks(**off):
    return merge_settings({"checks": {k: False for k in off}})


def test_turning_off_a_check_passes_regions_on(settings):
    erv = region("erv", [ladder(0.0)], taxid=39068)
    host = region("h", [ladder(0.9)])
    assert by_id(decide([erv, host], checks(retrovirus=1)))["erv"].action == "remove"
    d = by_id(decide([erv, host], checks(related_genome=1)))
    assert (d["h"].rule, d["h"].action) == ("6", "remove") and d["erv"].rule == "2"


def test_lineage_repeat_and_repeat_family_switches(settings):
    rep = RepeatEvidence(0.95, 0.95, "ACGTTGCA" * 4, lineage_copies=2, distant_copies=0)
    fam = RepeatEvidence(0.9, 0.9, "ACGTTGCA" * 4, lineage_copies=0, distant_copies=0)
    regions = [region("r", [ladder(0.0)], repeats=rep), region("f", [ladder(0.0)], repeats=fam)]
    d = by_id(decide(regions, settings))
    assert (d["r"].rule, d["f"].rule) == ("4", "5")
    d = by_id(decide(regions, checks(repeat_family=1)))
    assert (d["r"].rule, d["f"].rule) == ("4", "6")
    d = by_id(decide(regions, checks(lineage_repeat=1)))
    assert (d["r"].rule, d["f"].rule) == ("6", "6")


def test_all_alignment_checks_off(settings):
    from boardercollie.rules import TEST_OFF
    s = checks(related_genome=1, lineage_repeat=1, repeat_family=1)
    bad = [Control(f"c{i}", "host", [ladder(0.0)]) for i in range(10)]
    result = decide([region("h", [ladder(0.9)]), region("erv", taxid=39068)], s, controls=bad)
    assert result.status == "decided" and result.relative_alignment_test == TEST_OFF
    d = by_id(result)
    assert (d["h"].rule, d["h"].action) == ("6", "remove")
    assert (d["erv"].rule, d["erv"].action) == ("2", "keep")

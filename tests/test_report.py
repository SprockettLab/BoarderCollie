import csv
import json

from boardercollie.config import merge_settings
from boardercollie.decisions_io import write_decisions
from boardercollie.report import RULE_MEANING, write_host_report
from boardercollie.rules import Alignment, Control, Region, decide_genome


def region(rid, frac, length=1000):
    return Region(rid, "GCA_T.1", f"seq_{rid}", 1, length, "EXCLUDE", "prok:CFB group bacteria", "Bacteroides", 816,
                  "action_report", [Alignment("GCA_L.1", "ladder", frac, 0.9)])


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def test_host_report_rows_and_summary(tmp_path):
    s = merge_settings({})
    controls = [Control(f"{r}_C{i}", "host", [Alignment("GCA_L.1", "ladder", 0.95)], r) for r in ("a", "b") for i in range(3)]
    controls.append(Control("FOREIGN_X", "foreign", [Alignment("GCA_L.1", "ladder", 0.0)]))
    result = decide_genome("GCA_T.1", [region("a", 0.9), region("b", 0.0), region("c", 0.0, 250)], controls, s,
                           lambda t: set())
    dec = tmp_path / "decisions.tsv"
    write_decisions(result, dec, "hash1")
    run = tmp_path / "run.json"
    run.write_text(json.dumps({"removed_bp": 1250, "original_bp": 10000, "clean_bp": 8750, "original_md5": "m1",
                               "clean_md5": "m2", "run_key": "k1"}))
    rec = tmp_path / "record.json"
    rec.write_text(json.dumps({"tax_id": 1, "md5": "m1", "report": {"organism": {"organism_name": "Testus primus"}}}))
    out_r, out_g = tmp_path / "regions.tsv", tmp_path / "genomes.tsv"
    write_host_report("Test_host", [{"accession": "GCA_T.1", "decisions": str(dec), "run": str(run), "record": str(rec),
                                     "relative": "GCA_R.1", "ladder": ["GCA_L.1"], "screen": "run:2023-01-24"}],
                      out_r, out_g)
    r = {x["region_id"]: x for x in rows(out_r)}
    assert (r["a"]["action"], r["a"]["reason"]) == ("keep", RULE_MEANING["3"])
    assert (r["b"]["action"], r["b"]["reason"]) == ("remove", RULE_MEANING["6"])
    assert (r["c"]["action"], r["c"]["rule"]) == ("remove", "too_short")
    assert all(x["organism"] == "Testus primus" and x["host"] == "Test_host" for x in r.values())
    g = rows(out_g)[0]
    assert (g["regions_flagged"], g["regions_kept"], g["regions_removed"]) == ("3", "1", "2")
    assert (g["status"], g["relative_alignment_test"], g["host_controls"]) == ("decided", "applied", "6 of 6")
    assert (g["removed_bp"], g["relative"], g["ladder"], g["run_key"]) == ("1250", "GCA_R.1", "GCA_L.1", "k1")
    assert g["fcs_gx"] == "FCS-GX run by BoarderCollie (database build 2023-01-24)"


def test_download_only_genome(tmp_path):
    rec = tmp_path / "record.json"
    rec.write_text(json.dumps({"tax_id": 9606, "md5": "abc", "source": "local"}))
    out_r, out_g = tmp_path / "regions.tsv", tmp_path / "genomes.tsv"
    write_host_report("H", [{"accession": "x", "decisions": None, "run": str(rec), "record": str(rec)}], out_r, out_g)
    assert rows(out_r) == []
    g = rows(out_g)[0]
    assert g["status"] == "not screened (download_only)" and g["original_md5"] == "abc" and g["organism"] == "tax-id 9606"
    assert g["fcs_gx"] == "not screened"

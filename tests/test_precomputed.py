"""FCS-GX reports made earlier: finding them, checking them against the genome, and the workflow using them."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from boardercollie.config import ConfigError, parse_config
from boardercollie.fcs_report import ReportError
from boardercollie.fasta import write_fasta
from boardercollie.precomputed import build_date, check_reports, find_reports, header_info, reports_hash

from test_fcs_report import TAX_HEADER, tax_row

# Header details as FCS-GX 0.5.5 writes them (from a real report, 2026-11-06).
META = ('{"git-rev": "v0.5.5-9-gfc836c1", "run-date": "Fri Nov  6 12:32:20 2026", '
        '"db": {"build-date": "2023-01-24", "seqs": 3025824, "Gbp": 709.264}, '
        '"run-info": {"agg-cvg": 0.980024, "asserted-div": "anml:primates", '
        '"inferred-primary-divs": ["anml:primates", "anml:mammals", "anml:rodents"], '
        '"corrected-primary-divs": ["anml:primates", "anml:mammals", "anml:rodents"], '
        '"genome-hash": "608bfb626fac23ae"}}')
HEAD = '##[["FCS genome report", 2, 1], ' + META + ']\n'
COLS = "#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n"
TAX_HEAD = '##[["GX taxonomy analysis report", 3, 1], ' + META + ']\n'


def genome(tmp_path, name="G.1"):
    path = tmp_path / f"{name}.fna"
    write_fasta([("chr1", "ACGT" * 2500), ("bact1", "GGCC" * 750)], path, width=60)  # 10,000 and 3,000 bp
    return path


def reports(folder, key="G.1", action_rows=None, tax_rows=None, head=HEAD, tax_head=TAX_HEAD, taxid=""):
    folder.mkdir(parents=True, exist_ok=True)
    stem = f"{key}.{taxid}" if taxid else key
    a, t = folder / f"{stem}.fcs_gx_report.txt", folder / f"{stem}.taxonomy.rpt"
    a.write_text(head + COLS + "".join(action_rows if action_rows is not None
                                       else ["bact1\t1\t3000\t3000\tEXCLUDE\tprok:x\t100\tBacteroides\n"]))
    t.write_text(tax_head + TAX_HEADER + "\n" + "\n".join(tax_rows if tax_rows is not None else [
        tax_row("bact1", 3000, "Bacteroides", 816, "prok:x", "contaminant(div)"),
        tax_row("chr1~~1..5000", 10000, "Macaca fascicularis", 9541, "anml:primates", "primary-div"),
        tax_row("chr1~~5001..10000", 10000, "Macaca fascicularis", 9541, "anml:primates", "primary-div")]) + "\n")
    return a, t


# Finding the reports.

def test_find_both_namings(tmp_path):
    a, t = reports(tmp_path / "r")
    assert find_reports(tmp_path / "r", "G.1") == (a, t)
    a2, t2 = reports(tmp_path / "s", taxid="9541")  # as run_gx.py names them
    assert find_reports(tmp_path / "s", "G.1") == (a2, t2)


def test_find_none_and_no_prefix_match(tmp_path):
    reports(tmp_path / "r", key="G.10")
    assert find_reports(tmp_path / "r", "G.1") is None
    assert find_reports(tmp_path / "missing", "G.1") is None


def test_find_half_a_pair_or_two_pairs_is_an_error(tmp_path):
    a, t = reports(tmp_path / "r")
    t.unlink()
    with pytest.raises(ReportError, match="no matching taxonomy report"):
        find_reports(tmp_path / "r", "G.1")
    reports(tmp_path / "s")
    reports(tmp_path / "s", taxid="9541")
    with pytest.raises(ReportError, match="several"):
        find_reports(tmp_path / "s", "G.1")


# Header details.

def test_header_fields(tmp_path):
    a, t = reports(tmp_path / "r")
    assert header_info(a) == {"git-rev": "v0.5.5-9-gfc836c1", "build-date": "2023-01-24",
                              "asserted-div": "anml:primates", "genome-hash": "608bfb626fac23ae"}
    assert build_date(a, t) == "2023-01-24"
    a.write_text(a.read_text().replace("2023-01-24", "2024-05-01"))
    with pytest.raises(ReportError, match="different database builds"):
        build_date(a, t)


def test_header_without_details(tmp_path):
    a, t = reports(tmp_path / "r", head='##[["FCS genome report", 2, 1], {}]\n', tax_head="")
    assert header_info(a) == {} and header_info(t) == {}
    assert build_date(a, t) is None


def test_hash_changes_with_contents(tmp_path):
    a, t = reports(tmp_path / "r")
    h = reports_hash(a, t)
    t.write_text(t.read_text() + "\n")
    assert reports_hash(a, t) != h


# Checking against the genome.

def test_matching_reports_pass(tmp_path):
    a, t = reports(tmp_path / "r")
    rec = check_reports(a, t, genome(tmp_path), 9541)
    assert (rec["action_rows"], rec["taxonomy_rows"], rec["regions_to_test"], rec["genome_sequences"]) == (1, 3, 1, 2)
    assert rec["database_build"] == "2023-01-24" and rec["asserted_division"] == "anml:primates"
    assert rec["fcs_gx_genome_hash"] == "608bfb626fac23ae" and rec["fcs_gx_version"] == "v0.5.5-9-gfc836c1"


def test_sequence_missing_from_genome(tmp_path):
    a, t = reports(tmp_path / "r", action_rows=["chrX\t1\t3000\t3000\tEXCLUDE\tprok:x\t100\tBacteroides\n"])
    with pytest.raises(ReportError, match="chrX, which is not in the genome"):
        check_reports(a, t, genome(tmp_path), 9541)


def test_length_differs(tmp_path):
    a, t = reports(tmp_path / "r", action_rows=["bact1\t1\t3000\t3001\tEXCLUDE\tprok:x\t100\tBacteroides\n"])
    with pytest.raises(ReportError, match="length 3001"):
        check_reports(a, t, genome(tmp_path), 9541)
    a, t = reports(tmp_path / "s", tax_rows=[tax_row("chr1", 9999, "M", 9541, "anml:primates", "primary-div")])
    with pytest.raises(ReportError, match="length 9999"):
        check_reports(a, t, genome(tmp_path), 9541)
    a, t = reports(tmp_path / "u", tax_rows=[tax_row("chr1~~5001..10001", 10000, "M", 9541, "anml:primates", "x")])
    with pytest.raises(ReportError, match="ends past"):
        check_reports(a, t, genome(tmp_path), 9541)


def test_reports_from_two_runs(tmp_path):
    a, t = reports(tmp_path / "r")
    t.write_text(t.read_text().replace("608bfb626fac23ae", "0000000000000000"))
    with pytest.raises(ReportError, match="different FCS-GX runs"):
        check_reports(a, t, genome(tmp_path), 9541)


def test_headers_without_details_still_checked_by_sequence(tmp_path):
    a, t = reports(tmp_path / "s", head='##[["FCS genome report", 2, 1], {}]\n', tax_head="")
    rec = check_reports(a, t, genome(tmp_path), 9606)
    assert rec["asserted_division"] == "not stated" and rec["tax_id_checked"].startswith("no")


def test_empty_reports_pass(tmp_path):
    a, t = reports(tmp_path / "r", action_rows=[], tax_rows=[])
    rec = check_reports(a, t, genome(tmp_path), 9541)
    assert (rec["action_rows"], rec["regions_to_test"]) == (0, 0)


def test_malformed_report(tmp_path):
    a, t = reports(tmp_path / "r")
    a.write_text(a.read_text().replace("top_tax_name", "other"))
    with pytest.raises(ReportError, match="missing columns"):
        check_reports(a, t, genome(tmp_path), 9541)


# Config.

def test_config_reports_without_database():
    base = {"library_dir": "/l", "index_dir": "/i", "work_dir": "/w", "hosts": [{"name": "h", "genomes": ["G.1"]}]}
    c = parse_config({**base, "fcs_gx_reports": "/reports"})
    assert c.fcs_gx_reports == Path("/reports") and c.fcs_gx_db is None
    with pytest.raises(ConfigError, match="fcs_gx_db"):
        parse_config(base)


# The workflow (FCS-GX not installed: only the precomputed path can run).

SNAKEFILE = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"


def snakemake(cfg, cwd, *extra):
    pytest.importorskip("snakemake")
    env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}
    return subprocess.run([sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE), "--configfile", str(cfg),
                           "--cores", "1", *extra], cwd=cwd, env=env, capture_output=True, text=True)


def workflow_config(tmp_path, folder, **extra):
    fasta = genome(tmp_path, "host1")
    cfg = tmp_path / "config.yaml"
    off = {"related_genome": False, "lineage_repeat": False, "repeat_family": False}  # no network, no relatives
    cfg.write_text(yaml.safe_dump({
        "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
        "fcs_gx_reports": str(folder),
        "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 9541}], "indexes": []}],
        "settings": {"checks": off}, **extra}))
    return cfg


def test_workflow_uses_reports_and_never_screens(tmp_path):
    a, t = reports(tmp_path / "reports", key="host1")
    cfg = workflow_config(tmp_path, tmp_path / "reports")
    plan = snakemake(cfg, tmp_path, "--dry-run")
    assert plan.returncode == 0, plan.stderr
    out = plan.stdout + plan.stderr
    assert "rule screen_precomputed:" in out and "rule screen:" not in out
    assert "used for 1 of 1 genomes" in out
    done = snakemake(cfg, tmp_path, "--until", "screen_precomputed")
    assert done.returncode == 0, done.stderr
    d = tmp_path / "lib" / "genomes" / "host1" / "fcs_gx" / "precomputed" / reports_hash(a, t)
    assert (d / "fcs_gx_report.txt").read_text() == a.read_text()
    check = json.loads((d / "check.json").read_text())
    assert check["source"] == "precomputed" and check["genome_tax_id"] == 9541 and check["regions_to_test"] == 1


def test_workflow_stops_on_mismatched_reports(tmp_path):
    reports(tmp_path / "reports", key="host1",
            action_rows=["bact1\t1\t3000\t2999\tEXCLUDE\tprok:x\t100\tBacteroides\n"])
    done = snakemake(workflow_config(tmp_path, tmp_path / "reports"), tmp_path, "--until", "screen_precomputed")
    assert done.returncode != 0 and "length 2999" in done.stdout + done.stderr


def test_workflow_needs_database_for_genomes_without_reports(tmp_path):
    (tmp_path / "reports").mkdir()
    plan = snakemake(workflow_config(tmp_path, tmp_path / "reports"), tmp_path, "--dry-run")
    assert plan.returncode != 0 and "has no reports for: host1" in plan.stdout + plan.stderr


def test_workflow_screens_genomes_without_reports_where_it_always_did(tmp_path):
    db = tmp_path / "gxdb"
    db.mkdir()
    (db / "all.meta.jsonl").write_text('{"build-date": "2023-01-24"}\n')
    (tmp_path / "reports").mkdir()
    plan = snakemake(workflow_config(tmp_path, tmp_path / "reports", fcs_gx_db=str(db)), tmp_path, "--dry-run")
    out = plan.stdout + plan.stderr
    assert plan.returncode == 0, out
    assert "rule screen:" in out and "rule screen_precomputed:" not in out
    assert str(tmp_path / "lib" / "genomes" / "host1" / "fcs_gx" / "2023-01-24" / "fcs_gx_report.txt") in out


def test_workflow_needs_a_build_date(tmp_path):
    reports(tmp_path / "reports", key="host1", head='##[["FCS genome report", 2, 1], {}]\n', tax_head="")
    plan = snakemake(workflow_config(tmp_path, tmp_path / "reports"), tmp_path, "--dry-run")
    assert plan.returncode != 0 and "fcs_gx_db_build" in plan.stdout + plan.stderr
    ok = snakemake(workflow_config(tmp_path, tmp_path / "reports", fcs_gx_db_build="2023-01-24"), tmp_path, "--dry-run")
    assert ok.returncode == 0, ok.stderr

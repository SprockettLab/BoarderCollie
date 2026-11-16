"""A finished run is not redone by the next run with the same config.

Runs the workflow with download_only on a small local FASTA (no network, no FCS-GX),
then asks Snakemake, under a new run id, what it would do. It must be nothing:
each run has its own id, and a step whose command text changed would be rerun.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

pytest.importorskip("snakemake")

SNAKEFILE = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"


def snakemake(cfg, cwd, run_id, *extra):
    # Steps call "python -m boardercollie.cli": the interpreter running the tests must come first.
    env = {**os.environ, "BOARDERCOLLIE_RUN_ID": run_id,
           "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}
    return subprocess.run([sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE), "--configfile", str(cfg),
                           "--cores", "1", *extra], cwd=cwd, env=env, capture_output=True, text=True)


def test_second_run_redoes_nothing(tmp_path):
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTAC\n>c2\nGGGCCCAAAT\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(yaml.safe_dump({
        "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
        "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 9606}], "indexes": []}],
        "settings": {"download_only": True}}))
    first = snakemake(cfg, tmp_path, "first")
    assert first.returncode == 0, first.stderr
    assert (tmp_path / "lib" / "hosts" / "H").exists()
    second = snakemake(cfg, tmp_path, "second", "--dry-run")
    assert second.returncode == 0, second.stderr
    assert "Nothing to be done" in second.stdout + second.stderr, second.stdout + second.stderr
    log = (tmp_path / "work" / "logs" / "run_log.tsv").read_text()
    assert "\tfirst\tstart\tfetch\t" in log and "\tfirst\tfinish\tfetch\t" in log


def test_no_run_id_in_step_commands():
    text = SNAKEFILE.read_text()
    assert '"{RUN_ID}"' not in text and "{RUN_ID}" not in text.split("def logged")[1].split("def tool_versions")[0]


def test_each_genome_aligned_to_its_own_ladder(tmp_path):
    """A host joining genomes of two lineages: each genome's flagged regions are planned
    for alignment to its own relative and ladder only, never to the other genome's."""
    fa, fb = tmp_path / "primate.fna", tmp_path / "parasite.fna"
    fa.write_text(">a1\nACGTACGTAC\n")
    fb.write_text(">b1\nGGGCCCAAAT\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(yaml.safe_dump({
        "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
        "fcs_gx_db": str(tmp_path / "gxdb"), "fcs_gx_db_build": "2023-01-24",
        "hosts": [{"name": "H", "indexes": [], "genomes": [
            {"fasta": str(fa), "tax_id": 9544, "relative": "GCA_100.1", "ladder": ["GCF_101.1", "GCF_102.1"]},
            {"fasta": str(fb), "tax_id": 5850, "relative": "GCA_200.1", "ladder": ["GCF_201.1", "GCF_202.1"]}]}]}))
    dry = snakemake(cfg, tmp_path, "dry", "--dry-run", "--printshellcmds")
    out = dry.stdout + dry.stderr
    assert dry.returncode == 0, out
    planned = {}
    for line in out.splitlines():
        for acc in ("primate", "parasite"):
            for ref in ("GCA_100.1", "GCF_101.1", "GCF_102.1", "GCA_200.1", "GCF_201.1", "GCF_202.1"):
                if f"/{acc}/" in line and f"{ref}.paf" in line:
                    planned.setdefault(acc, set()).add(ref)
    assert planned.get("primate") and planned.get("parasite"), out[-3000:]
    assert planned["primate"] <= {"GCA_100.1", "GCF_101.1", "GCF_102.1"}, planned
    assert planned["parasite"] <= {"GCA_200.1", "GCF_201.1", "GCF_202.1"}, planned
    assert "GCF_102.1" in planned["primate"] and "GCF_202.1" in planned["parasite"]  # each distant reference


def test_adding_a_parasite_leaves_the_primate_genome_alone(tmp_path):
    """The primate genome's run key and every file planned for it are the same whether
    the host is the primate alone or the primate joined with a parasite: a cleaned
    primate genome is reused, not cleaned again, when the parasite is added."""
    import re
    fa, fb = tmp_path / "primate.fna", tmp_path / "parasite.fna"
    fa.write_text(">a1\nACGTACGTAC\n")
    fb.write_text(">b1\nGGGCCCAAAT\n")
    primate = {"fasta": str(fa), "tax_id": 9544, "relative": "GCA_100.1", "ladder": ["GCF_101.1", "GCF_102.1"]}
    parasite = {"fasta": str(fb), "tax_id": 5850, "relative": "GCA_200.1", "ladder": ["GCF_201.1", "GCF_202.1"]}

    def planned(name, genomes):
        d = tmp_path / name
        d.mkdir()
        cfg = d / "config.yaml"
        cfg.write_text(yaml.safe_dump({
            "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"),
            "work_dir": str(tmp_path / "work"), "fcs_gx_db": str(tmp_path / "gxdb"),
            "fcs_gx_db_build": "2023-01-24", "hosts": [{"name": "H", "indexes": [], "genomes": genomes}]}))
        dry = snakemake(cfg, d, name, "--dry-run", "--printshellcmds")
        out = dry.stdout + dry.stderr
        assert dry.returncode == 0, out
        files = set(re.findall(r"[^\s\"';,]*/primate/[^\s\"';,]*", out))
        keys = set(re.findall(r"/genomes/primate/runs/([0-9a-f]{12})/", out))
        return files, keys

    alone_files, alone_keys = planned("alone", [primate])
    joined_files, joined_keys = planned("joined", [primate, parasite])
    assert len(alone_keys) == 1 and alone_keys == joined_keys
    assert alone_files and alone_files == joined_files

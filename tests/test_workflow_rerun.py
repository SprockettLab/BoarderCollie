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

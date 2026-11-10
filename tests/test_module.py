"""BoarderCollie loaded as a Snakemake module by another workflow, as MAGmaker does.

The calling workflow names a fixed path for the host reference (reference_out) and
asks only for that file. Runs with download_only on a small local FASTA (no network,
no FCS-GX).
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from boardercollie.config import ConfigError, parse_config

pytest.importorskip("snakemake")

REPO = Path(__file__).resolve().parent.parent

PARENT = """
module boardercollie:
    snakefile: "{snakefile}"
    config: config["boardercollie"]

use rule * from boardercollie as boardercollie_*

rule all:
    default_target: True
    input: config["host_filter"]["genome"]
"""


def run(cwd, *extra):
    env = {**os.environ, "BOARDERCOLLIE_RUN_ID": "test",
           "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}
    return subprocess.run([sys.executable, "-m", "snakemake", "--configfile", "config.yaml", "--cores", "1", *extra],
                          cwd=cwd, env=env, capture_output=True, text=True)


def test_module_writes_host_reference_to_named_path(tmp_path):
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTAC\n>c2\nGGGCCCAAAT\n")
    (tmp_path / "Snakefile").write_text(PARENT.format(snakefile=REPO / "workflow" / "Snakefile"))
    genome = tmp_path / "host" / "H.fna"
    (tmp_path / "config.yaml").write_text(yaml.safe_dump({
        "host_filter": {"genome": str(genome)},
        "boardercollie": {
            "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
            "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 9606}], "indexes": [],
                       "reference_out": str(genome)}],
            "settings": {"download_only": True}}}))
    first = run(tmp_path)
    assert first.returncode == 0, first.stdout + first.stderr
    assert genome.read_text() == fasta.read_text()
    members = (tmp_path / "host" / "H.fna.members.tsv").read_text().splitlines()
    assert members[0].startswith("accession\t") and members[1].startswith("host1\t")
    # Resources (foreign controls) and the package come from BoarderCollie's folder, not the caller's.
    assert not (tmp_path / "resources").exists()
    again = run(tmp_path, "--dry-run")
    assert again.returncode == 0 and "Nothing to be done" in again.stdout + again.stderr, again.stdout + again.stderr


def test_reference_out_config():
    base = {"library_dir": "/l", "index_dir": "/i", "work_dir": "/w", "fcs_gx_db": "/db"}
    c = parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "reference_out": "/x/h.fna"}]})
    assert c.hosts[0].reference_out == "/x/h.fna"
    assert parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"]}]}).hosts[0].reference_out is None
    with pytest.raises(ConfigError, match="file path"):
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "reference_out": ""}]})
    with pytest.raises(ConfigError, match="another host"):
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "reference_out": "/x.fna"},
                                         {"name": "k", "genomes": ["G.2"], "reference_out": "/x.fna"}]})
    with pytest.raises(ConfigError, match="unknown keys"):
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "refrence_out": "/x.fna"}]})

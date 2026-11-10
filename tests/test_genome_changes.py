"""A host whose genomes change gets a new reference, and the old one is kept.

Runs the workflow with download_only on small local FASTA files (no network, no FCS-GX).
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from boardercollie.config import write_if_changed

pytest.importorskip("snakemake")

SNAKEFILE = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"
GENOME_1 = ">c1\nACGTACGTAC\n>c2\nGGGCCCAAAT\n"
GENOME_2 = ">d1\nTTGACCATGG\n>d2\nCAGTTGACCA\n"


def run(tmp_path, genomes, *extra):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(yaml.safe_dump({
        "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
        "hosts": [{"name": "H", "genomes": [{"fasta": str(g), "tax_id": 9606} for g in genomes], "indexes": []}],
        "settings": {"download_only": True}}))
    env = {**os.environ, "BOARDERCOLLIE_RUN_ID": "test",
           "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}
    result = subprocess.run([sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE), "--configfile", str(cfg),
                             "--cores", "1", *extra], cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def host_dirs(tmp_path):
    return {p.name: (p / "H.fna").read_text() for p in (tmp_path / "lib" / "hosts" / "H").iterdir()}


def test_other_genome_same_settings_gets_its_own_host_key(tmp_path):
    a, b = tmp_path / "a" / "host.fna", tmp_path / "b" / "other.fna"
    for path, text in ((a, GENOME_1), (b, GENOME_2)):
        path.parent.mkdir()
        path.write_text(text)
    run(tmp_path, [a])
    first = host_dirs(tmp_path)
    assert list(first.values()) == [GENOME_1]
    # Same settings, other genome: before host keys named the genomes, this wrote
    # over the first reference in the same folder.
    run(tmp_path, [b])
    second = host_dirs(tmp_path)
    assert len(second) == 2
    assert sorted(second.values()) == sorted([GENOME_1, GENOME_2])
    assert {k: second[k] for k in first} == first


def test_local_fasta_replaced_at_same_path_is_copied_again(tmp_path):
    fasta = tmp_path / "host.fna"
    fasta.write_text(GENOME_1)
    run(tmp_path, [fasta])
    library_copy = tmp_path / "lib" / "genomes" / "host" / "original" / "host.fna"
    assert library_copy.read_text() == GENOME_1
    # The library copy is made older than the replaced file, as it is when a file is
    # replaced later; moving the source into the future instead looks like clock skew.
    old = library_copy.stat().st_mtime - 60
    os.utime(library_copy, (old, old))
    fasta.write_text(GENOME_2)
    run(tmp_path, [fasta])
    assert library_copy.read_text() == GENOME_2
    assert GENOME_2 in host_dirs(tmp_path).values()
    again = run(tmp_path, [fasta], "--dry-run")
    assert "Nothing to be done" in again.stdout + again.stderr


def test_config_used_is_not_rewritten_when_unchanged(tmp_path):
    fasta = tmp_path / "host.fna"
    fasta.write_text(GENOME_1)
    run(tmp_path, [fasta])
    used = tmp_path / "work" / "config.used.yaml"
    before = used.stat()
    run(tmp_path, [fasta], "--dry-run")
    after = used.stat()
    assert (before.st_ino, before.st_mtime_ns) == (after.st_ino, after.st_mtime_ns)


def test_write_if_changed(tmp_path):
    path = tmp_path / "sub" / "f.yaml"
    assert write_if_changed(path, "a: 1\n")
    assert path.read_text() == "a: 1\n"
    inode = path.stat().st_ino
    assert not write_if_changed(path, "a: 1\n")
    assert path.stat().st_ino == inode
    assert write_if_changed(path, "a: 2\n")
    assert path.read_text() == "a: 2\n"
    assert [p.name for p in path.parent.iterdir()] == ["f.yaml"]
    assert path.stat().st_mode & 0o044

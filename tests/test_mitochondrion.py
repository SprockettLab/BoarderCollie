"""A host genome's mitochondrion: the config key, the checks, the join, the workflow plan.

The checks run against stand-ins for the NCBI lookups (no network).
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from boardercollie.cli import cmd_combine
from boardercollie.config import ConfigError, parse_config
from boardercollie.mitochondrion import MitochondrionError, fetch_mitochondrion

REPO = Path(__file__).resolve().parent.parent
BASE = {"library_dir": "/l", "index_dir": "/i", "work_dir": "/w", "fcs_gx_db": "/db"}
SEQ = "ACGT" * 4000 + "A"  # 16,001 bp

# Alouatta caraya (9502) and A. seniculus (9503) are species; 99999 is a made-up subspecies of 9502.
RANKS = {9502: {"species": 9502, "genus": 9499}, 9503: {"species": 9503, "genus": 9499},
         99999: {"species": 9502, "genus": 9499}}


def summary(acc, **over):
    s = {"accessionversion": acc, "title": "Alouatta caraya mitochondrion, complete genome", "taxid": 9502,
         "organism": "Alouatta caraya", "slen": len(SEQ), "genome": "mitochondrion", "completeness": "complete",
         "sourcedb": "refseq"}
    s.update(over)
    return s


def run(tmp_path, acc="NC_064185.1", genome="GCA_059577405.1", record=None, held=(), seq=SEQ, ncbi=None):
    record = record or {"source": "ncbi", "tax_id": 9502}
    asked = []

    def assembly_mitochondria(g):
        asked.append(g)
        return list(held)
    rec = fetch_mitochondrion(acc, genome, record, tmp_path / "mt.fna",
                              summary=lambda a: summary(a, **(ncbi or {})),
                              fasta=lambda a: f">{a} Alouatta caraya mitochondrion\n{seq}\n",
                              assembly_mitochondria=assembly_mitochondria,
                              ranks_of=lambda ids: {i: RANKS[i] for i in ids})
    return rec, asked


def test_config_key():
    c = parse_config({**BASE, "hosts": [{"name": "h", "genomes": [
        {"accession": "GCA_059577405.1", "mitochondrion": "NC_064185.1"}, "GCA_000001405.29"]}]})
    assert c.hosts[0].genomes[0].mitochondrion == "NC_064185.1"
    assert c.hosts[0].genomes[1].mitochondrion is None
    for bad in ("NC_064185", "auto", "", 5):
        with pytest.raises(ConfigError, match="mitochondrion must be"):
            parse_config({**BASE, "hosts": [{"name": "h", "genomes": [
                {"accession": "GCA_059577405.1", "mitochondrion": bad}]}]})


def test_written_with_its_record(tmp_path):
    rec, asked = run(tmp_path)
    text = (tmp_path / "mt.fna").read_text().splitlines()
    assert text[0].startswith(">NC_064185.1 ") and "".join(text[1:]) == SEQ
    assert asked == ["GCA_059577405.1"]
    assert rec["length"] == len(SEQ) and rec["species_tax_id"] == 9502 and rec["genome"] == "GCA_059577405.1"
    assert "holds no mitochondrion" in rec["assembly_check"]


def test_subspecies_record_of_the_same_species(tmp_path):
    rec, _ = run(tmp_path, ncbi={"taxid": 99999})
    assert rec["tax_id"] == 99999 and rec["species_tax_id"] == 9502


def test_other_species_refused(tmp_path):
    with pytest.raises(MitochondrionError, match="not from the species"):
        run(tmp_path, ncbi={"taxid": 9503, "organism": "Alouatta seniculus"})
    assert not (tmp_path / "mt.fna").exists()


@pytest.mark.parametrize("over", [{"genome": "genomic"}, {"completeness": "partial"}])
def test_not_a_complete_mitochondrion_refused(tmp_path, over):
    with pytest.raises(MitochondrionError, match="not a complete mitochondrial genome"):
        run(tmp_path, ncbi=over)


def test_assembly_that_already_holds_one_refused(tmp_path):
    with pytest.raises(MitochondrionError, match="already holds a mitochondrion"):
        run(tmp_path, held=["CM045128.1"])
    assert not (tmp_path / "mt.fna").exists()


def test_length_must_match_ncbi(tmp_path):
    with pytest.raises(MitochondrionError, match="NCBI states"):
        run(tmp_path, seq=SEQ[:-1])


def test_local_fasta_genome_not_asked(tmp_path):
    rec, asked = run(tmp_path, record={"source": "local", "tax_id": 9502})
    assert asked == [] and rec["assembly_check"].startswith("not checked")


def test_combine_lists_what_each_member_is(tmp_path):
    g, m = tmp_path / "g.fna", tmp_path / "m.fna"
    g.write_text(">chr1\nACGT\n")
    m.write_text(">NC_064185.1 mito\nGGCC\n")

    class Args:
        accessions, fastas, records = ["G.1", "NC_064185.1"], [str(g), str(m)], ["rg", "rm"]
        member_of = ["-", "G.1"]
        out, members_out = str(tmp_path / "h.fna"), str(tmp_path / "members.tsv")
    cmd_combine(Args)
    assert (tmp_path / "h.fna").read_text() == ">chr1\nACGT\n>NC_064185.1 mito\nGGCC\n"
    rows = [l.split("\t") for l in (tmp_path / "members.tsv").read_text().splitlines()]
    assert rows[0][-1] == "member_of" and rows[1][0] == "G.1" and rows[1][-1] == "-"
    assert rows[2][0] == "NC_064185.1" and rows[2][-1] == "G.1"


def plan(tmp_path, mitochondrion):
    pytest.importorskip("snakemake")
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTAC\n")
    genome = {"fasta": str(fasta), "tax_id": 9502}
    if mitochondrion:
        genome["mitochondrion"] = mitochondrion
    cfg = tmp_path / "config.yaml"
    cfg.write_text(yaml.safe_dump({
        "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
        "hosts": [{"name": "H", "genomes": [genome], "indexes": []}], "settings": {"download_only": True}}))
    env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}
    p = subprocess.run([sys.executable, "-m", "snakemake", "-s", str(REPO / "workflow" / "Snakefile"),
                        "--configfile", str(cfg), "--cores", "1", "--dry-run", "--printshellcmds"],
                       cwd=tmp_path, env=env, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr


def host_key_of(out):
    keys = {l.split("/hosts/H/")[1].split("/")[0] for l in out.splitlines() if "/hosts/H/" in l}
    assert len(keys) == 1, keys
    return keys.pop()


def test_workflow_plans_and_joins_the_mitochondrion(tmp_path):
    with_mt = plan(tmp_path, "NC_064185.1")
    assert "rule fetch_mitochondrion:" in with_mt
    assert "fetch-mitochondrion --accession NC_064185.1 --genome host1 " in with_mt
    assert "--accessions host1 NC_064185.1 " in with_mt and "--member-of - host1 " in with_mt


def test_host_key_changes_only_with_a_mitochondrion(tmp_path):
    a, b, c = (tmp_path / x for x in "abc")
    for d in (a, b, c):
        d.mkdir()
    without, again, with_mt = plan(a, None), plan(b, None), plan(c, "NC_064185.1")
    assert "rule fetch_mitochondrion:" not in without and "--member-of - " in without
    # A host key names genomes and settings, not folders, so plans from different folders compare.
    assert host_key_of(without) == host_key_of(again)
    assert host_key_of(with_mt) != host_key_of(without)


def test_assembly_shared_by_two_hosts_only_one_with_a_mitochondrion(tmp_path):
    # One assembly can serve two hosts (a species with none of its own borrows a sister
    # species'); only the host that names a mitochondrion gets it.
    pytest.importorskip("snakemake")
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTAC\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(yaml.safe_dump({
        "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
        "hosts": [{"name": "A", "genomes": [{"fasta": str(fasta), "tax_id": 9502, "mitochondrion": "NC_064185.1"}],
                   "indexes": []},
                  {"name": "B", "genomes": [{"fasta": str(fasta), "tax_id": 9502}], "indexes": []}],
        "settings": {"download_only": True}}))
    env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")}
    p = subprocess.run([sys.executable, "-m", "snakemake", "-s", str(REPO / "workflow" / "Snakefile"),
                        "--configfile", str(cfg), "--cores", "1", "--dry-run", "--printshellcmds"],
                       cwd=tmp_path, env=env, capture_output=True, text=True)
    out = p.stdout + p.stderr
    assert p.returncode == 0, out
    combines = [l for l in out.splitlines() if "cli combine " in l]
    assert any("--accessions host1 NC_064185.1 " in l and "/hosts/A/" in l for l in combines), combines
    assert any("--accessions host1 --fastas" in l and "/hosts/B/" in l for l in combines), combines

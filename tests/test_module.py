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


def test_index_out_config():
    base = {"library_dir": "/l", "index_dir": "/i", "work_dir": "/w", "fcs_gx_db": "/db"}
    c = parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "indexes": ["bowtie2", "minimap2"],
                                         "bowtie2_out": "/x/h", "minimap2_out": "/x/h.mmi"}]})
    with pytest.raises(ConfigError, match="needs bowtie2 in indexes"):   # indexes are off by default
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "bowtie2_out": "/x/h"}]})
    assert (c.hosts[0].bowtie2_out, c.hosts[0].minimap2_out) == ("/x/h", "/x/h.mmi")
    assert parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"]}]}).hosts[0].bowtie2_out is None
    with pytest.raises(ConfigError, match="needs bowtie2 in indexes"):
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "indexes": ["minimap2"],
                                         "bowtie2_out": "/x/h"}]})
    with pytest.raises(ConfigError, match="needs minimap2 in indexes"):
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "indexes": [],
                                         "minimap2_out": "/x/h.mmi"}]})
    with pytest.raises(ConfigError, match="another host or output"):
        parse_config({**base, "hosts": [{"name": "h", "genomes": ["G.1"], "reference_out": "/x/h",
                                         "bowtie2_out": "/x/h"}]})


INDEX_PARENT = """
module boardercollie:
    snakefile: "{snakefile}"
    config: config["boardercollie"]

use rule * from boardercollie as boardercollie_*

rule all:
    default_target: True
    input: config["want"]
"""


def test_index_out_built_once_and_shared(tmp_path):
    # Two calling workflows in separate folders, same host: the indexes are built once in
    # index_dir; each folder gets hard links at its own fixed paths. The second builds nothing.
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\n" + "ACGTTGCAAGGCTTAACCGGTATCGATCGGCTAGCTAGGATCCA" * 30 + "\n")
    shared = {"library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"),
              "work_dir": str(tmp_path / "work")}
    outs = []
    for run_dir in ("runA", "runB", "runC"):
        d = tmp_path / run_dir
        d.mkdir()
        (d / "Snakefile").write_text(INDEX_PARENT.format(snakefile=REPO / "workflow" / "Snakefile"))
        stem, mmi = d / "out" / "H", d / "out" / "H.map-ont.mmi"
        (d / "config.yaml").write_text(yaml.safe_dump({
            "want": [str(d / "out" / "H.fna"), str(stem) + ".index_done.txt", str(mmi)],
            "boardercollie": {**shared, "hosts": [{
                "name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 9606}],
                "indexes": ["bowtie2", "minimap2"], "reference_out": str(d / "out" / "H.fna"),
                "bowtie2_out": str(stem), "minimap2_out": str(mmi)}],
                "settings": {"download_only": True}}}))
        # No checksums: Snakemake skips the mtime check for files under 1 MB whose
        # checksum is unchanged, which a real host reference (gigabytes) never gets.
        done = run(d, "--max-checksum-file-size", "0")
        assert done.returncode == 0, done.stdout + done.stderr
        outs.append(done.stdout + done.stderr)
        built = sorted((tmp_path / "idx" / "bt2" / "H").glob("*/H.1.bt2l"))
        assert len(built) == 1, built
        for p in ("1", "2", "3", "4", "rev.1", "rev.2"):
            link = Path(f"{stem}.{p}.bt2l")
            assert link.stat().st_ino == (built[0].parent / f"H.{p}.bt2l").stat().st_ino, p
        assert mmi.stat().st_ino == next((tmp_path / "idx" / "minimap2" / "H").glob("*/H.map-ont.mmi")).stat().st_ino
        assert Path(f"{stem}.index_done.txt").read_text().startswith("shared_index\t")
    assert "rule boardercollie_bowtie2:" in outs[0] and "rule boardercollie_minimap2_index:" in outs[0]
    # Each caller's reference_out (a hard link of the library FASTA) is touched by Snakemake
    # after its job, which bumps the library file's mtime; later callers must still reuse.
    for i, name in ((1, "second"), (2, "third")):
        assert "rule boardercollie_bowtie2:" not in outs[i], f"{name} workflow rebuilt the bowtie2 index"
        assert "rule boardercollie_minimap2_index:" not in outs[i], f"{name} workflow rebuilt the minimap2 index"
        assert "rule boardercollie_bowtie2_out_H:" in outs[i] and "rule boardercollie_minimap2_out_H:" in outs[i]
    # And the first caller, re-planned after the others ran, has nothing to do.
    again = run(tmp_path / "runA", "--dry-run", "--max-checksum-file-size", "0")
    assert again.returncode == 0 and "Nothing to be done" in again.stdout + again.stderr, again.stdout + again.stderr


# Step commands are written out when the plan is made, so a dry run with commands printed (-p)
# fails on any command naming something only BoarderCollie's own Snakefile defines. Run as a
# module, those names are not in reach: a command must take them from its params.
# download_only skips every such step, so these two plans run with the checks on.

def module_plan(tmp_path, bc):
    (tmp_path / "Snakefile").write_text(PARENT.format(snakefile=REPO / "workflow" / "Snakefile"))
    genome = tmp_path / "host" / "H.fna"
    bc = {"library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"),
          "work_dir": str(tmp_path / "work"), **bc}
    bc["hosts"][0]["reference_out"] = str(genome)
    (tmp_path / "config.yaml").write_text(yaml.safe_dump({"host_filter": {"genome": str(genome)},
                                                          "boardercollie": bc}))
    plan = run(tmp_path, "--dry-run", "--printshellcmds")
    return plan.returncode, plan.stdout + plan.stderr


def test_module_plan_with_every_check(tmp_path):
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTACGGGCCCAAAT\n")
    gxdb = tmp_path / "gxdb"
    gxdb.mkdir()
    (gxdb / "all.meta.jsonl").write_text('{"build-date": "2023-01-24"}\n')
    # Relative and ladder named, so the plan needs no NCBI lookup.
    code, out = module_plan(tmp_path, {
        "fcs_gx_db": str(gxdb),
        "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 60711, "relative": "GCF_000409795.2"}],
                   "ladder": ["GCF_049350105.2", "GCF_029281585.2"], "indexes": []}]})
    assert code == 0, out
    for step in ("screen", "extract", "flagged", "queries", "align", "copies_own", "copies_ref", "trf",
                 "decide", "clean", "combine", "reference_out_H"):
        assert f"rule boardercollie_{step}:" in out, step
    work = tmp_path / "work"
    assert f"extract --config {work / 'config.used.yaml'} " in out
    assert f"decide --config {work / 'config.used.yaml'} " in out
    assert out.count("--settings-hash ") == 2 and "--settings-hash {" not in out


def _every_check_plan(tmp_path, controls_fna_text):
    """Plan with every check on, a shared controls FASTA older than the controls table (as
    for any caller whose BoarderCollie clone is new), and the given FASTA contents."""
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTACGGGCCCAAAT\n")
    gxdb = tmp_path / "gxdb"
    gxdb.mkdir()
    (gxdb / "all.meta.jsonl").write_text('{"build-date": "2023-01-24"}\n')
    work = tmp_path / "work"
    work.mkdir()
    fna = work / "foreign_controls.fna"
    fna.write_text(controls_fna_text)
    os.utime(fna, (1_000_000_000, 1_000_000_000))   # 2001: older than the clone's table
    return module_plan(tmp_path, {
        "fcs_gx_db": str(gxdb),
        "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 60711, "relative": "GCF_000409795.2"}],
                   "ladder": ["GCF_049350105.2", "GCF_029281585.2"], "indexes": []}]})


def test_fresh_clone_does_not_refetch_current_controls(tmp_path):
    from boardercollie.controls import listed
    accs = listed(REPO / "resources" / "foreign_controls.tsv")
    code, out = _every_check_plan(tmp_path, "".join(f">FOREIGN_{a}\nACGT\n" for a in accs))
    assert code == 0, out
    assert "rule boardercollie_foreign_controls:" not in out, out


def test_changed_controls_list_is_fetched_again(tmp_path):
    code, out = _every_check_plan(tmp_path, ">FOREIGN_NOT_LISTED.1\nACGT\n")
    assert code == 0, out
    assert "rule boardercollie_foreign_controls:" in out, out


def test_controls_up_to_date(tmp_path):
    from boardercollie.controls import fetched, listed, up_to_date
    tsv = tmp_path / "c.tsv"
    tsv.write_text("# accession\tdescription\nA.1\tone\nB.2\ttwo\n")
    fna = tmp_path / "c.fna"
    assert listed(tsv) == ["A.1", "B.2"] and not up_to_date(tsv, fna)
    fna.write_text(">FOREIGN_A.1\nACGT\n>FOREIGN_B.2\nGGCC\n")
    assert fetched(fna) == ["A.1", "B.2"] and up_to_date(tsv, fna)
    fna.write_text(">FOREIGN_A.1\nACGT\n")
    assert not up_to_date(tsv, fna)


def test_module_plan_with_reports_made_earlier(tmp_path):
    from test_precomputed import genome as small_genome, reports
    fasta = small_genome(tmp_path, "host1")
    reports(tmp_path / "reports", key="host1")
    code, out = module_plan(tmp_path, {
        "fcs_gx_reports": str(tmp_path / "reports"),
        "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 9541}], "indexes": []}],
        "settings": {"checks": {"related_genome": False, "lineage_repeat": False, "repeat_family": False}}})
    assert code == 0, out
    assert "rule boardercollie_screen_precomputed:" in out and " --build 2023-01-24 " in out


def test_steps_find_the_package_without_pythonpath(tmp_path):
    # With --use-conda, Snakemake runs every step with PYTHONPATH removed, so a step's
    # "python -m boardercollie.cli" found no package (ModuleNotFoundError on demon). The
    # caller here removes PYTHONPATH before each step the same way; the steps must still run.
    fasta = tmp_path / "host1.fna"
    fasta.write_text(">c1\nACGTACGTAC\n>c2\nGGGCCCAAAT\n")
    (tmp_path / "Snakefile").write_text('shell.prefix("set -euo pipefail; unset PYTHONPATH; ")\n'
                                        + PARENT.format(snakefile=REPO / "workflow" / "Snakefile"))
    genome = tmp_path / "host" / "H.fna"
    (tmp_path / "config.yaml").write_text(yaml.safe_dump({
        "host_filter": {"genome": str(genome)},
        "boardercollie": {
            "library_dir": str(tmp_path / "lib"), "index_dir": str(tmp_path / "idx"), "work_dir": str(tmp_path / "work"),
            "hosts": [{"name": "H", "genomes": [{"fasta": str(fasta), "tax_id": 9606}], "indexes": [],
                       "reference_out": str(genome)}],
            "settings": {"download_only": True}}}))
    done = run(tmp_path)
    assert done.returncode == 0, done.stdout + done.stderr
    assert genome.read_text() == fasta.read_text()

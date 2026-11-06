"""The workflow's Python steps on small synthetic genomes and alignments."""

import random

import pytest

from boardercollie.clean import check_accounting, filtered_report
from boardercollie.cli import cmd_combine, lineage_members
from boardercollie.composition import cpg_oe, gc_fraction, percentile
from boardercollie.config import merge_settings
from boardercollie.decide import decide_from_files, read_tsv
from boardercollie.extract import extract
from boardercollie.fasta import IndexedFasta, read_fasta, write_fasta
from boardercollie.paf import count_copies, read_paf, summarize
from boardercollie.trf import read_trf_ngs

from test_fcs_report import TAX_HEADER, tax_row

ACTION_HEAD = ('##[["FCS genome report", 2, 1], {}]\n'
               "#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n")


def rand_seq(n, rng, gc=0.4):
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


@pytest.fixture
def genome(tmp_path):
    rng = random.Random(1)
    seqs = [("chr1", rand_seq(20000, rng)), ("chr2", rand_seq(15000, rng)),
            ("bact1", rand_seq(3000, rng, 0.6)), ("chr3", rand_seq(9000, rng))]
    # chr3 is half N: windows there must be redrawn.
    seqs[3] = ("chr3", "N" * 6000 + seqs[3][1][6000:])
    path = tmp_path / "genome.fna"
    write_fasta(seqs, path, width=60)
    action = tmp_path / "report.txt"
    action.write_text(ACTION_HEAD + "bact1\t1\t3000\t3000\tEXCLUDE\tprok:x\t100\tBacteroides\n"
                      "chr2\t5001\t6000\t15000\tFIX\tprok:x\t7\tPrevotella\n")
    tax = tmp_path / "tax.rpt"
    tax.write_text("\n".join([TAX_HEADER,
                              tax_row("bact1", 3000, "Bacteroides", 816, "prok:x", "contaminant(div)"),
                              tax_row("chr2~~5001..6000", 15000, "Prevotella", 838, "prok:x", "contaminant(div)"),
                              tax_row("chr1", 20000, "Macaca fascicularis", 9541, "anml:primates", "primary-div")]) + "\n")
    return path, action, tax, dict(seqs)


def test_indexed_fasta(genome, tmp_path):
    path, _, _, seqs = genome
    fa = IndexedFasta(path, tmp_path / "g.fai")
    assert fa.lengths()["chr1"] == 20000
    assert fa.fetch("chr2", 59, 125) == seqs["chr2"][58:125]
    assert fa.fetch("bact1", 1, 3000) == seqs["bact1"]
    with pytest.raises(ValueError):
        fa.fetch("chr1", 0, 10)


def test_composition():
    assert gc_fraction("GGCCAATT") == 0.5
    assert abs(cpg_oe("CGCGCGCG") - 2.0) < 1e-9
    assert cpg_oe("AAAA") == 0.0
    assert percentile([1, 2, 3, 4, 5], 50) == 3
    assert percentile([1, 2, 3, 4, 5], 0) == 1 and percentile([1, 2, 3, 4, 5], 100) == 5


def test_extract(genome, tmp_path):
    path, action, tax, seqs = genome
    s = merge_settings({"background_per_segment": 20})
    out = [tmp_path / n for n in ("seg.fna", "seg.tsv", "comp.tsv")]
    n = extract(path, action, tax, s, tmp_path / "g.fai", *out)
    assert n == 2
    rows = read_tsv(out[1])
    flagged = [r for r in rows if r["type"] == "flagged"]
    controls = [r for r in rows if r["type"] == "host_control"]
    assert len(flagged) == 2 and len(controls) == 6
    # Controls never come from sequences with a flagged region, and have the flagged length.
    assert {r["seq_id"] for r in controls} <= {"chr1", "chr3"}
    assert all(float(r["n_frac"]) <= 0.1 for r in controls)
    assert sorted(int(r["length"]) for r in controls) == [1000] * 3 + [3000] * 3
    fasta = read_fasta(out[0])
    assert fasta["R2"] == seqs["chr2"][5000:6000]
    # Same seed, same controls.
    out2 = [tmp_path / n for n in ("s2.fna", "s2.tsv", "c2.tsv")]
    extract(path, action, tax, s, tmp_path / "g.fai", *out2)
    assert out[1].read_text() == out2[1].read_text()
    assert len(read_tsv(out[2])) == 2


def test_paf_summary_and_copies(tmp_path):
    paf = tmp_path / "a.paf"
    paf.write_text("\n".join([
        "q1\t1000\t0\t400\t+\tt1\t50000\t20000\t20400\t380\t400\t60",
        "q1\t1000\t300\t600\t+\tt2\t700\t100\t400\t250\t300\t60",
        "q2\t500\t0\t500\t+\tt3\t600\t50\t550\t500\t500\t60",
        "q2\t500\t0\t450\t+\tsrc\t9000\t999\t1449\t450\t450\t0",
    ]) + "\n")
    s = summarize(read_paf(paf), {"q1": 1000, "q2": 500, "q3": 100}, flank_bp=10000)
    assert abs(s["q1"].aligned_fraction - 0.6) < 1e-9 and s["q1"].embedded
    assert s["q2"].aligned_fraction == 1.0 and not s["q2"].embedded
    assert s["q3"].aligned_fraction == 0.0
    c = count_copies(read_paf(paf), 0.5, {"q2": ("src", 1000, 1499)})
    assert c["q1"] == 0 and c["q2"] == 1


def test_trf(tmp_path):
    p = tmp_path / "t.ngs"
    p.write_text("@R1\n1 600 16 37.5 16 95 0 900 25 25 25 25 2.0 ACGTTGCAAGGCTTAC ACGT\n"
                 "@R2\n")
    t = read_trf_ngs(p, {"R1": 1000, "R2": 500, "R3": 50})
    assert (t["R1"].cover, t["R1"].main_period, t["R1"].main_unit) == (0.6, 16, "ACGTTGCAAGGCTTAC")
    assert t["R2"].cover == 0 and t["R3"].main_unit == ""


def test_decide_from_files_end_to_end(genome, tmp_path):
    path, action, tax, _ = genome
    s = merge_settings({"background_per_segment": 10})
    seg, segt, comp = tmp_path / "seg.fna", tmp_path / "seg.tsv", tmp_path / "comp.tsv"
    extract(path, action, tax, s, tmp_path / "g.fai", seg, segt, comp)
    foreign = tmp_path / "foreign.fna"
    write_fasta([("FOREIGN_X", "ACGT" * 500)], foreign)
    controls = [r for r in read_tsv(segt) if r["type"] == "host_control"]
    # Host controls align fully to the ladder; R1 (bact1) aligns nowhere; R2 (chr2 piece) aligns to the ladder.
    lines = [f"{r['id']}\t{r['length']}\t0\t{r['length']}\t+\tL\t100000\t500\t{500 + int(r['length'])}\t{r['length']}\t{r['length']}\t60"
             for r in controls] + ["R2\t1000\t0\t900\t+\tL\t100000\t500\t1400\t880\t900\t60"]
    ladder = tmp_path / "ladder.paf"
    ladder.write_text("\n".join(lines) + "\n")
    rel = tmp_path / "rel.paf"
    rel.write_text("")
    empty = tmp_path / "empty.paf"
    empty.write_text("")
    trf = tmp_path / "trf.ngs"
    trf.write_text("@R1\n@R2\n")
    out, rem = tmp_path / "decisions.tsv", tmp_path / "removals.tsv"
    res = decide_from_files("GCA_TEST.1", segt, seg, comp, foreign,
                            [("REL.1", "relative", rel), ("L", "ladder", ladder)], trf, empty, [], empty,
                            s, tmp_path / "tax.json", out, rem, "hash",
                            fetch_taxonomy=lambda ids: {i: {"parents": set(), "ranks": {}, "name": ""} for i in ids})
    assert res.status == "decided"
    d = {x.region.seq_id: x for x in res.decisions}
    assert (d["bact1"].rule, d["bact1"].action) == ("6", "remove")
    assert (d["chr2"].rule, d["chr2"].action) == ("3", "keep")
    removals = read_tsv(rem)
    assert [(r["seq_id"], r["start"], r["end"]) for r in removals] == [("bact1", "1", "3000")]


def test_filtered_report_and_accounting(genome, tmp_path):
    path, action, _, seqs = genome
    rem = tmp_path / "removals.tsv"
    rem.write_text("region_id\tseq_id\tstart\tend\tlength\tseq_len\tfcs_action\tsource\n"
                   "R1\tbact1\t1\t3000\t3000\t3000\tEXCLUDE\taction_report\n"
                   "R3\tchr1\t101\t200\t100\t20000\tnone\ttaxonomy_report\n")
    out = tmp_path / "applied.txt"
    assert filtered_report(action, rem, out) == 2
    rows = [l.split("\t") for l in out.read_text().splitlines() if not l.startswith("#")]
    assert [r[0] for r in rows] == ["bact1", "chr1"] and rows[1][4] == "FIX"
    # Accounting: simulate gx clean-genome output.
    clean, removed = tmp_path / "clean.fna", tmp_path / "removed.fna"
    write_fasta([("chr1_a", seqs["chr1"][:100]), ("chr1_b", seqs["chr1"][200:]), ("chr2", seqs["chr2"]),
                 ("chr3", seqs["chr3"])], clean)
    write_fasta([("bact1", seqs["bact1"]), ("chr1_x", seqs["chr1"][100:200])], removed)
    acc = check_accounting(path, clean, removed, rem)
    assert acc["removed_bp"] == 3100 and acc["clean_bp"] + 3100 == acc["original_bp"]
    write_fasta([("bact1", seqs["bact1"])], removed)
    with pytest.raises(ValueError):
        check_accounting(path, clean, removed, rem)


def test_filtered_report_missing_row(genome, tmp_path):
    _, action, _, _ = genome
    rem = tmp_path / "removals.tsv"
    rem.write_text("region_id\tseq_id\tstart\tend\tlength\tseq_len\tfcs_action\tsource\n"
                   "R9\tchrX\t1\t10\t10\t10\tEXCLUDE\taction_report\n")
    with pytest.raises(ValueError):
        filtered_report(action, rem, tmp_path / "o.txt")


def test_lineage_members(tmp_path):
    def rec(name, taxid):
        p = tmp_path / f"{name}.json"
        p.write_text(f'{{"tax_id": {taxid}}}')
        return p
    own = rec("own", 9534)
    cands = [f"A={rec('a', 60711)}=a.paf", f"B={rec('b', 9606)}=b.paf"]
    tax = {9534: {"ranks": {"genus": 9533}}, 60711: {"ranks": {"genus": 9533}}, 9606: {"ranks": {"genus": 9605}}}
    keep = lineage_members(own, cands, "genus", fetch=lambda ids: tax)
    assert [a for a, _ in keep] == ["A"]


def test_combine_rejects_duplicate_names(tmp_path):
    a, b = tmp_path / "a.fna", tmp_path / "b.fna"
    write_fasta([("s1", "ACGT")], a)
    write_fasta([("s1", "ACGT")], b)

    class Args:
        accessions, fastas, records = ["A", "B"], [str(a), str(b)], ["ra", "rb"]
        out, members_out = str(tmp_path / "h.fna"), str(tmp_path / "m.tsv")
    with pytest.raises(ValueError):
        cmd_combine(Args)


def test_taxonomy_cache_written_with_nothing_to_look_up(tmp_path):
    from boardercollie.decide import lineage_lookup
    cache = tmp_path / "tax.json"
    lookup = lineage_lookup([], cache, fetch=lambda ids: {})
    assert cache.exists() and lookup(816) == set()

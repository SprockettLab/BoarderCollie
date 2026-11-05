"""Reading FCS-GX reports, with small reports in FCS-GX 0.5.5's format."""

import pytest

from boardercollie.fcs_report import ReportError, read_candidates

ACTION = """##[["FCS genome report", 2, 1], {"db": {"build-date": "2023-01-24"}}]
#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name
ctgA\t1\t5000\t5000\tEXCLUDE\tprok:a-proteobacteria\t100\tPhyllobacterium calauticae
ctgB\t2001\t3000\t9000\tFIX\tprok:CFB group bacteria\t11\tPrevotella sp.
"""

TAX_HEADER = "#seq-id\tseq-len\t(xp,lc,co,n,mt,pt,pm)-len\tcvg-by-all\tsep1\ttax-name-1\ttax-id-1\tdiv-1\t" + \
    "\t".join(f"x{i}" for i in range(9, 31)) + "\tresult\tdiv\tdiv_pct_cvg"


def tax_row(seq_id, seq_len, name, taxid, div, result):
    return "\t".join([seq_id, str(seq_len), "0,0,0,0,0,0,0", "0", "|", name, str(taxid), div] +
                     ["."] * 22 + [result, div, "100"])


TAX = "\n".join([
    '##[["GX taxonomy analysis report", 3, 1], {}]',
    TAX_HEADER,
    tax_row("ctgA", 5000, "Phyllobacterium calauticae", 2817027, "prok:a-proteobacteria", "contaminant(div)"),
    tax_row("ctgB~~1..2000", 9000, "Macaca fascicularis", 9541, "anml:primates", "primary-div"),
    tax_row("ctgB~~2001..3000", 9000, "Prevotella sp.", 838, "prok:CFB group bacteria", "contaminant(div)"),
    tax_row("ctgC~~98167..99926", 120000, "Phyllobacterium calauticae", 2817027, "prok:a-proteobacteria", "contaminant(div)"),
    tax_row("ctgD", 4000, "Macaca fascicularis", 9541, "anml:primates", "primary-div"),
]) + "\n"


@pytest.fixture
def reports(tmp_path):
    a, t = tmp_path / "g.fcs_gx_report.txt", tmp_path / "g.taxonomy.rpt"
    a.write_text(ACTION)
    t.write_text(TAX)
    return a, t


def test_action_rows_get_taxids(reports):
    c = {x.seq_id: x for x in read_candidates(*reports)}
    assert (c["ctgA"].fcs_action, c["ctgA"].fcs_top_taxid, c["ctgA"].length) == ("EXCLUDE", 2817027, 5000)
    assert (c["ctgB"].start, c["ctgB"].end, c["ctgB"].fcs_top_taxid) == (2001, 3000, 838)


def test_taxonomy_only_contaminant_piece_is_a_candidate(reports):
    c = {x.seq_id: x for x in read_candidates(*reports)}
    x = c["ctgC"]
    assert (x.source, x.fcs_action, x.start, x.end) == ("taxonomy_report", "none", 98167, 99926)


def test_host_pieces_are_not_candidates(reports):
    assert "ctgD" not in {x.seq_id for x in read_candidates(*reports)}
    assert len(read_candidates(*reports)) == 3


def test_missing_column_is_an_error(tmp_path, reports):
    bad = tmp_path / "bad.txt"
    bad.write_text(ACTION.replace("top_tax_name", "something_else"))
    with pytest.raises(ReportError):
        read_candidates(bad, reports[1])

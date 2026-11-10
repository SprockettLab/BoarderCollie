from boardercollie.config import merge_settings
from boardercollie.decisions_io import COLUMNS, write_decisions
from boardercollie.rules import Alignment, Control, Region, decide_genome


def test_decisions_file(tmp_path):
    r = Region("r1", "GCA_TEST.1", "ctgA", 1, 500, "EXCLUDE", "prok:x", "Bacteroides", 816, "action_report",
               [Alignment("GCA_LADDER.1", "ladder", 0.0)])
    controls = [Control("c", "host", [Alignment("GCA_LADDER.1", "ladder", 0.9)])]
    result = decide_genome("GCA_TEST.1", [r], controls, merge_settings({}), lambda t: set())
    out = tmp_path / "decisions.tsv"
    write_decisions(result, out, "abc123")
    lines = out.read_text().splitlines()
    comments = [l for l in lines if l.startswith("#")]
    assert "# relative-genome alignment test: applied" in comments
    assert "# settings hash: abc123" in comments
    header, row = [l.split("\t") for l in lines if not l.startswith("#")]
    assert tuple(header) == COLUMNS and len(row) == len(COLUMNS)
    assert row[COLUMNS.index("rule")] == "6" and row[COLUMNS.index("action")] == "remove"

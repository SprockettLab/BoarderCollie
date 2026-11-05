from pathlib import Path

import pytest

from boardercollie.config import DEFAULT_SETTINGS, ConfigError, load_config, merge_settings, parse_config

EXAMPLE = Path(__file__).resolve().parents[1] / "config" / "config.example.yaml"


def minimal(**extra):
    d = {"library_dir": "/lib", "index_dir": "/idx", "work_dir": "/work", "fcs_gx_db": "/db",
         "hosts": [{"name": "Chlorocebus_spp", "genomes": ["GCF_047675955.1", "GCA_023783515.1"]}]}
    d.update(extra)
    return d


def test_example_config_loads():
    c = load_config(EXAMPLE)
    assert c.hosts and c.settings["host_like_min"] == DEFAULT_SETTINGS["host_like_min"]


def test_defaults_applied():
    c = parse_config(minimal())
    assert c.settings == merge_settings({})
    assert c.hosts[0].indexes == ["bowtie2", "minimap2"]


def test_override_and_hash_changes():
    a = parse_config(minimal())
    b = parse_config(minimal(settings={"repeat_trf_min": 0.7}))
    assert b.settings["repeat_trf_min"] == 0.7
    assert a.settings_hash() != b.settings_hash()
    assert a.settings_hash() == parse_config(minimal()).settings_hash()


@pytest.mark.parametrize("settings", [
    {"no_such_setting": 1},
    {"host_like_min": 1.5},
    {"foreign_like_max": 0.6},
    {"on_controls_failed": "maybe"},
    {"repeat_lineage_rank": "kingdom"},
    {"keep_taxa": ["Retroviridae"]},
    {"conflict_background_percentiles": [99, 1]},
])
def test_bad_settings_rejected(settings):
    with pytest.raises(ConfigError):
        parse_config(minimal(settings=settings))


def test_local_fasta_needs_taxid():
    with pytest.raises(ConfigError):
        parse_config(minimal(hosts=[{"name": "x", "genomes": [{"fasta": "/g.fna"}]}]))
    c = parse_config(minimal(hosts=[{"name": "x", "genomes": [{"fasta": "/g.fna", "tax_id": 9534}]}]))
    assert c.hosts[0].genomes[0].tax_id == 9534


def test_missing_and_duplicate():
    d = minimal()
    del d["fcs_gx_db"]
    with pytest.raises(ConfigError):
        parse_config(d)
    with pytest.raises(ConfigError):
        parse_config(minimal(hosts=[{"name": "x", "genomes": ["A.1", "A.1"]}]))

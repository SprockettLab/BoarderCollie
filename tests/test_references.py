import json

from boardercollie.references import Assembly, cached_choice, choose_relative, is_automatic

RANKS = {"species": 1, "genus": 10, "family": 100, "order": 1000, "class": 10000}


def rep(acc, taxid, project="P1", sample=None, level="Scaffold", ref=False, n50=1000, paired=None):
    return {"accession": acc, "paired_accession": paired, "organism": {"tax_id": taxid, "organism_name": f"sp{taxid}"},
            "assembly_info": {"bioproject_accession": project, "biosample": {"accession": sample or f"S_{acc}"},
                              "assembly_level": level, "refseq_category": "reference genome" if ref else None},
            "assembly_stats": {"contig_n50": n50, "scaffold_n50": n50}}


SCREENED = Assembly.from_report(rep("GCA_1.1", 1, project="P1", sample="S1", paired="GCF_1.1"))


def lookup(table):
    calls = []

    def f(taxon):
        calls.append(taxon)
        return table.get(taxon, [])
    f.calls = calls
    return f


def test_automatic_markers():
    for v in (None, "", "none", "NA", "na", " None ", "auto"):
        assert is_automatic(v)
    assert not is_automatic("GCA_000001405.29")


def test_self_twin_and_same_animal_excluded():
    f = lookup({1: [rep("GCA_1.1", 1), rep("GCF_1.1", 1, paired="GCA_1.1"), rep("GCA_2.1", 1, sample="S1")],
                10: [rep("GCA_3.1", 2)]})
    ch = choose_relative(SCREENED, RANKS, f)
    assert ch["relative"] == "GCA_3.1" and ch["rank"] == "genus"
    why = {r["accession"]: r["excluded_because"] for r in ch["candidates"] if not r["eligible"]}
    assert set(why) == {"GCA_1.1", "GCA_2.1"} or set(why) == {"GCF_1.1", "GCA_2.1"}
    assert "BioSample" in why["GCA_2.1"]


def test_deepest_rank_wins_over_quality():
    f = lookup({1: [rep("GCA_5.1", 1, project="P9", level="Contig", n50=10)],
                10: [rep("GCF_6.1", 2, project="P8", level="Chromosome", ref=True, n50=10**8)]})
    ch = choose_relative(SCREENED, RANKS, f)
    assert ch["relative"] == "GCA_5.1" and f.calls == [1]


def test_different_project_preferred_at_same_rank():
    f = lookup({10: [rep("GCA_7.1", 2, project="P1", level="Chromosome", ref=True, n50=10**8),
                     rep("GCA_8.1", 3, project="P2", level="Contig", n50=10)]})
    ch = choose_relative(SCREENED, RANKS, f)
    assert ch["relative"] == "GCA_8.1"
    assert "different BioProject" in ch["reason"]


def test_same_project_used_when_only_choice():
    f = lookup({10: [rep("GCA_7.1", 2, project="P1")]})
    ch = choose_relative(SCREENED, RANKS, f)
    assert ch["relative"] == "GCA_7.1" and "same BioProject" in ch["reason"]


def test_quality_order_then_refseq_twin():
    f = lookup({10: [rep("GCA_9.1", 2, project="P2", level="Scaffold", ref=True, n50=5),
                     rep("GCA_10.1", 2, project="P2", level="Chromosome", n50=10**7),
                     rep("GCA_11.1", 2, project="P2", level="Chromosome", ref=True, n50=100, paired="GCF_11.1"),
                     rep("GCF_11.1", 2, project="P2", level="Chromosome", ref=True, n50=100, paired="GCA_11.1")]})
    ch = choose_relative(SCREENED, RANKS, f)
    assert ch["relative"] == "GCF_11.1"
    assert [r["accession"] for r in ch["candidates"]].count("GCA_11.1") == 0


def test_nothing_found():
    ch = choose_relative(SCREENED, RANKS, lookup({}))
    assert ch["relative"] is None and "no eligible assembly" in ch["reason"]


def test_choice_is_cached(tmp_path):
    p = tmp_path / "choice.json"
    first = cached_choice(p, lambda: {"relative": "GCA_3.1"})
    again = cached_choice(p, lambda: {"relative": "GCA_99.1"})
    assert first["relative"] == again["relative"] == "GCA_3.1"
    assert "chosen_on" in json.loads(p.read_text())


class FakeNCBI:
    """Stands in for boardercollie.ncbi: species 1 in genus 10, family 100, order 1000."""

    def __init__(self):
        self.calls = 0
        self.lineages = {1: RANKS, 2: {**RANKS, "species": 2},
                         3: {**RANKS, "species": 3, "genus": 11},
                         4: {**RANKS, "species": 4, "genus": 12, "family": 101}}
        self.reports = {"GCA_1.1": rep("GCA_1.1", 1, project="P1", sample="S1", n50=10**7)}
        self.under = {1: [rep("GCA_1.1", 1, project="P1", sample="S1", n50=10**7)],
                      10: [rep("GCA_2.1", 2, project="P2", n50=10**7), rep("GCA_20.1", 2, project="P3", n50=10)],
                      100: [rep("GCF_3.1", 3, project="P4", ref=True, level="Chromosome", n50=10**8)],
                      1000: [rep("GCF_4.1", 4, project="P5", ref=True, level="Chromosome", n50=10**8)]}

    def genome_report(self, acc):
        self.calls += 1
        return self.reports[acc]

    def ranks_of(self, taxids):
        self.calls += 1
        return {t: self.lineages[t] for t in taxids}

    def genomes_under(self, taxid, reference_only=False):
        self.calls += 1
        rows = [r for t, rs in self.under.items() for r in rs
                if self.lineages[r["organism"]["tax_id"]].get(
                    next(k for k, v in RANKS.items() if v == taxid), None) == taxid] if taxid in RANKS.values() else []
        if reference_only:
            rows = [r for r in rows if r["assembly_info"]["refseq_category"] == "reference genome"]
        return rows


def config(tmp_path, relative=None, ladder=None, settings=None):
    from boardercollie.config import parse_config
    h = {"name": "host", "genomes": [{"accession": "GCA_1.1", "relative": relative}]}
    if ladder is not None:
        h["ladder"] = ladder
    return parse_config({"library_dir": str(tmp_path), "index_dir": "/i", "work_dir": "/w", "fcs_gx_db": "/d",
                         "hosts": [h], "settings": settings or {}})


def test_resolve_fills_blanks_and_reuses_saved_choice(tmp_path):
    from boardercollie.references import resolve
    f = FakeNCBI()
    cfg = config(tmp_path)
    lines = resolve(cfg, tmp_path, f)
    g, h = cfg.hosts[0].genomes[0], cfg.hosts[0]
    # Genus-level relative with scaffold N50 >= 1 Mb (GCA_20.1 is too fragmented).
    assert g.relative == "GCA_2.1"
    assert h.ladder == ["GCF_3.1", "GCF_4.1"]  # family (other genus), order (other family)
    assert len(lines) == 2
    assert (tmp_path / "genomes" / "GCA_1.1" / "relative_choice.json").exists()
    assert (tmp_path / "hosts" / "host" / "ladder_choice.json").exists()
    calls = f.calls
    cfg2 = config(tmp_path)
    resolve(cfg2, tmp_path, f)
    assert f.calls == calls  # saved choices reused, NCBI not asked again
    assert cfg2.hosts[0].genomes[0].relative == "GCA_2.1" and cfg2.hosts[0].ladder == ["GCF_3.1", "GCF_4.1"]


def test_resolve_chooses_again_when_settings_change(tmp_path):
    from boardercollie.references import resolve
    f = FakeNCBI()
    resolve(config(tmp_path), tmp_path, f)
    cfg = config(tmp_path, settings={"relative_min_scaffold_n50": 0, "ladder_ranks": ["family"]})
    resolve(cfg, tmp_path, f)
    assert cfg.hosts[0].ladder == ["GCF_3.1"]


def test_resolve_keeps_named_references(tmp_path):
    from boardercollie.references import resolve
    f = FakeNCBI()
    cfg = config(tmp_path, relative="GCA_9.1", ladder=["GCA_8.1"])
    assert resolve(cfg, tmp_path, f) == [] and f.calls == 0
    assert cfg.hosts[0].genomes[0].relative == "GCA_9.1" and cfg.hosts[0].ladder == ["GCA_8.1"]


def test_resolve_records_none_found(tmp_path):
    from boardercollie.references import resolve
    f = FakeNCBI()
    f.under = {1: f.under[1]}
    cfg = config(tmp_path)
    lines = resolve(cfg, tmp_path, f)
    assert cfg.hosts[0].genomes[0].relative is None and cfg.hosts[0].ladder == []
    assert "none found" in lines[0] and "none found" in lines[1]


# ── Per-genome ladders: a host that joins genomes of different lineages ──────

PARASITE = {"species": 5, "genus": 50, "family": 500, "order": 5000, "class": 50000}


class MixedNCBI(FakeNCBI):
    """FakeNCBI plus a second lineage: species 5 in genus 50, family 500, order 5000,
    with a relative (species 6, same genus) and one reference genome per ladder rank."""

    def __init__(self):
        super().__init__()
        self.lineages.update({5: PARASITE, 6: {**PARASITE, "species": 6},
                              7: {**PARASITE, "species": 7, "genus": 51},
                              8: {**PARASITE, "species": 8, "genus": 52, "family": 501}})
        self.reports["GCA_5.1"] = rep("GCA_5.1", 5, project="Q1", sample="T1", n50=10**7)
        self.under.update({5: [self.reports["GCA_5.1"]],
                           50: [rep("GCA_6.1", 6, project="Q2", n50=10**7)],
                           500: [rep("GCF_7.1", 7, project="Q3", ref=True, level="Chromosome", n50=10**8)],
                           5000: [rep("GCF_8.1", 8, project="Q4", ref=True, level="Chromosome", n50=10**8)]})

    def genomes_under(self, taxid, reference_only=False):
        self.calls += 1
        rows = [r for rs in self.under.values() for r in rs
                if taxid in self.lineages[r["organism"]["tax_id"]].values()]
        rows = list({r["accession"]: r for r in rows}.values())
        if reference_only:
            rows = [r for r in rows if r["assembly_info"]["refseq_category"] == "reference genome"]
        return rows


def mixed_config(tmp_path, host_ladder=None, genome_ladder=None):
    from boardercollie.config import parse_config
    second = {"accession": "GCA_5.1"}
    if genome_ladder is not None:
        second["ladder"] = genome_ladder
    h = {"name": "host_and_parasite", "genomes": [{"accession": "GCA_1.1"}, second]}
    if host_ladder is not None:
        h["ladder"] = host_ladder
    return parse_config({"library_dir": str(tmp_path), "index_dir": "/i", "work_dir": "/w", "fcs_gx_db": "/d",
                         "hosts": [h]})


def test_mixed_host_gets_one_ladder_per_lineage(tmp_path):
    from boardercollie.references import resolve
    f = MixedNCBI()
    cfg = mixed_config(tmp_path)
    resolve(cfg, tmp_path, f)
    host, parasite = cfg.hosts[0].genomes
    # Each genome's ladder comes from its own lineage, never from the other's.
    assert host.ladder == ["GCF_3.1", "GCF_4.1"]
    assert parasite.ladder == ["GCF_7.1", "GCF_8.1"]
    assert host.relative == "GCA_2.1" and parasite.relative == "GCA_6.1"
    assert cfg.hosts[0].ladder == []  # no single host ladder
    d = tmp_path / "hosts" / "host_and_parasite"
    assert (d / "lineage_groups.json").exists()
    assert (d / "ladder_choice" / "GCA_1.1.json").exists() and (d / "ladder_choice" / "GCA_5.1.json").exists()
    assert not (d / "ladder_choice.json").exists()
    # Saved choices are reused: NCBI is not asked again.
    calls = f.calls
    cfg2 = mixed_config(tmp_path)
    resolve(cfg2, tmp_path, f)
    assert f.calls == calls
    assert [g.ladder for g in cfg2.hosts[0].genomes] == [["GCF_3.1", "GCF_4.1"], ["GCF_7.1", "GCF_8.1"]]


def test_single_lineage_host_unchanged(tmp_path):
    """A host whose genomes share a lineage gets the same ladder, in the same saved file
    with the same inputs, as before genomes had their own: run keys cannot change."""
    from boardercollie.references import resolve
    f = FakeNCBI()
    cfg = config(tmp_path)
    resolve(cfg, tmp_path, f)
    g, h = cfg.hosts[0].genomes[0], cfg.hosts[0]
    assert g.ladder == h.ladder == ["GCF_3.1", "GCF_4.1"]
    saved = json.loads((tmp_path / "hosts" / "host" / "ladder_choice.json").read_text())
    assert saved["inputs"] == {"genomes": ["GCA_1.1"], "relatives": ["GCA_2.1"],
                               "ladder_ranks": ["family", "order"], "relative_min_scaffold_n50": 1000000}
    assert not (tmp_path / "hosts" / "host" / "ladder_choice").exists()


def test_saved_single_ladder_used_without_ncbi(tmp_path):
    """A ladder saved by an earlier version (no lineage_groups.json beside it) is used
    as it is: the genomes' lineages are not looked up."""
    from boardercollie.references import resolve
    f = FakeNCBI()
    resolve(config(tmp_path), tmp_path, f)
    (tmp_path / "hosts" / "host" / "lineage_groups.json").unlink()
    calls = f.calls
    cfg = config(tmp_path)
    resolve(cfg, tmp_path, f)
    assert f.calls == calls
    assert cfg.hosts[0].genomes[0].ladder == ["GCF_3.1", "GCF_4.1"]


def test_genome_ladder_overrides_host_ladder(tmp_path):
    from boardercollie.references import resolve
    f = MixedNCBI()
    cfg = mixed_config(tmp_path, host_ladder=["GCA_90.1"], genome_ladder=["GCF_7.1", "GCF_8.1"])
    lines = resolve(cfg, tmp_path, f)
    assert lines and all(": relative " in l for l in lines)  # relatives chosen, no ladder chosen
    assert not (tmp_path / "hosts" / "host_and_parasite").exists()
    host, parasite = cfg.hosts[0].genomes
    assert host.ladder == ["GCA_90.1"] and parasite.ladder == ["GCF_7.1", "GCF_8.1"]


def test_named_ladder_for_one_genome_automatic_for_the_other(tmp_path):
    from boardercollie.references import resolve
    f = MixedNCBI()
    cfg = mixed_config(tmp_path, genome_ladder=["GCF_7.1", "GCF_8.1"])
    resolve(cfg, tmp_path, f)
    host, parasite = cfg.hosts[0].genomes
    assert parasite.ladder == ["GCF_7.1", "GCF_8.1"]
    assert host.ladder == ["GCF_3.1", "GCF_4.1"]


def test_missing_rank_keeps_one_group():
    """A rank one genome has no taxon for is left out, as choose_ladder leaves it out,
    so genomes that differ only there stay together."""
    from boardercollie.references import lineage_groups
    a = {"family": 100, "order": 1000}
    b = {"family": None, "order": 1000}
    c = {"family": 101, "order": 1000}
    assert lineage_groups(["a", "b", "c"], {"a": a, "b": b, "c": c}.get, ["family", "order"]) == [["a", "b", "c"]]
    d = {"family": 500, "order": 5000}
    assert lineage_groups(["a", "d"], {"a": a, "d": d}.get, ["family", "order"]) == [["a"], ["d"]]


def test_genome_ladder_config_checks():
    import pytest
    from boardercollie.config import ConfigError, parse_config

    def cfg(genome, ladder=None):
        h = {"name": "h", "genomes": [genome, "GCA_2.1"]}
        if ladder is not None:
            h["ladder"] = ladder
        return parse_config({"library_dir": "/l", "index_dir": "/i", "work_dir": "/w", "fcs_gx_db": "/d",
                             "hosts": [h]})

    c = cfg({"accession": "GCA_1.1", "ladder": ["GCF_9.1"]}, ladder=["GCF_8.1"])
    assert [g.ladder for g in c.hosts[0].genomes] == [["GCF_9.1"], ["GCF_8.1"]]
    assert cfg({"accession": "GCA_1.1", "ladder": "auto"}).hosts[0].genomes[0].ladder == []
    with pytest.raises(ConfigError):
        cfg({"accession": "GCA_1.1", "ladder": ["GCA_2.1"]})          # a ladder genome is a host genome
    with pytest.raises(ConfigError):
        cfg({"accession": "GCA_1.1", "relative": "GCF_9.1", "ladder": ["GCF_9.1"]})
    with pytest.raises(ConfigError):
        cfg({"accession": "GCA_1.1", "ladder": "GCF_9.1"})            # not a list

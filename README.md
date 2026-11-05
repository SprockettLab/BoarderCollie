# BoarderCollie

Host read removal in metagenomics aligns reads to a host genome and discards whatever aligns.
If the host genome assembly itself contains foreign DNA (bacteria, fungi, parasites, other animals, laboratory contaminants), real microbial reads are discarded as host.
BoarderCollie screens host genome assemblies for foreign DNA, tests every flagged region, removes only what fails, and builds the host references and indexes used for read removal.

Status: in development; not yet ready for use.

## How it works

1. **Screen** each genome with NCBI FCS-GX.
2. **Test** every flagged region with the relative-genome alignment test: does the region align to related host genomes as well as the same genome's own host DNA does? Random host windows (host controls) and known foreign sequences (foreign controls) calibrate the test for each genome.
3. **Decide** by written rules with configurable thresholds:
   - **C**: the genome's host controls must align and its foreign controls must not; otherwise the relative-genome alignment test is not applied to that genome, and the outputs say so.
   - **R0**: regions whose best FCS-GX match lies in Retroviridae (endogenous retroviruses) are kept.
   - **R1**: regions that align to related host genomes are kept.
   - **R2**: tandem repeats with copies in other genomes of the same lineage, and none in a distant reference, are kept.
   - **R2b**: tandem repeats whose repeat unit matches one kept by R1 or R2 are kept.
   - **R3**: everything else is removed.
4. **Clean** each genome with FCS-GX's own cleaning tool, applying only the removals.
5. **Combine** the cleaned genomes of a host into one reference and build bowtie2 and minimap2 indexes.

Every region, test result, rule and threshold is recorded, and the records can be written out as supplementary tables.

## Configuration

See [`config/config.example.yaml`](config/config.example.yaml).
A host is a name and a list of NCBI assembly accessions; every threshold is a setting with a default.

## Development

```
python -m venv .venv
.venv/bin/pip install -e ".[test]"
.venv/bin/pytest
```

## License

MIT; see [LICENSE](LICENSE).

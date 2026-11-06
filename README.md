# BoarderCollie

Host read removal in metagenomics aligns reads to a host genome and discards whatever aligns.
If the host genome assembly itself contains foreign DNA (bacteria, fungi, parasites, other animals, laboratory contaminants), real microbial reads are discarded as host.
BoarderCollie screens host genome assemblies for foreign DNA, tests every flagged region, removes only what fails, and builds the host references and indexes used for read removal.

Status: in development; not yet ready for use.

## How it works

1. **Screen** each genome with NCBI [FCS-GX](https://github.com/ncbi/fcs-gx).
2. **Test** every flagged region with the relative-genome alignment test: does the region align to related host genomes as well as the same genome's own host DNA does? Random host windows (host controls) and known foreign sequences (foreign controls) calibrate the test for each genome.
3. **Decide** which genomic regions are kept and which are discarded by the following configurable rules:
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
A host is a name, its genomes (NCBI assembly accessions, each with its closest relative for now) and a reference ladder; every threshold is a setting with a default.

## Running

```
mamba env create -f workflow/envs/boardercollie.yaml -p ./env
./env/bin/pip install --no-deps -e .
./env/bin/snakemake -s workflow/Snakefile --configfile config.yaml --cores 64 --resources mem_mb=740000
```

FCS-GX needs its database (about 520 GB) in memory, so screening runs one genome at a time on a large-memory node.
It is much faster when the database is on local storage: copy it once per run with `scripts/stage_fcs_gx_db.sh` and pass the copy with `--config fcs_gx_db=<copy>`, deleting the copy afterwards.

Outputs, under `library_dir`:

- `genomes/<accession>/original/`: the genome as downloaded and its NCBI record (md5 checked against NCBI).
- `genomes/<accession>/fcs_gx/<database build>/`: FCS-GX reports.
- `genomes/<accession>/runs/<run key>/`: `decisions.tsv`, the cleaned and removed sequence, the applied FCS-GX report and `run.json`.
- `hosts/<name>/<host key>/`: the joined host reference and `members.tsv`.

and under `index_dir`: `bt2/<name>/<host key>/` and `minimap2/<name>/<host key>/`.

## Development

```
python -m venv .venv
.venv/bin/pip install -e ".[test]"
.venv/bin/pytest
```

## License

MIT; see [LICENSE](LICENSE).

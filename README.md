<p align="center">
  <img src="images/boardercollie_logo.png" alt="BoarderCollie logo" width="240">
</p>

# BoarderCollie

Host read removal in metagenomics aligns reads to a host genome and discards whatever aligns.
If the host genome assembly itself contains foreign DNA (bacteria, fungi, parasites, other animals, laboratory contaminants), real microbial reads are discarded as host.
BoarderCollie screens host genome assemblies for foreign DNA, tests every flagged region, removes only what fails, and builds the host references and indexes used for read removal.

Status: in development; not yet ready for use.

## How it works

1. **Screen** each genome with NCBI [FCS-GX](https://github.com/ncbi/fcs-gx).
2. **Test** every flagged region with the relative-genome alignment test: does the region align to related host genomes as well as the same genome's own host DNA does? Random host windows (host controls) and known foreign sequences (foreign controls) calibrate the test for each genome.
3. **Decide** which genomic regions are kept and which are discarded by the following configurable rules:
   1. The genome's host controls must align and its foreign controls must not; otherwise the relative-genome alignment test is not applied to that genome, and the outputs say so.
   2. Regions whose best FCS-GX match lies in Retroviridae (endogenous retroviruses) are kept.
   3. Regions that align to related host genomes are kept.
   4. Tandem repeats with copies in other genomes of the same lineage, and none in a distant reference, are kept.
   5. Tandem repeats whose repeat unit matches one kept by rule 3 or 4 are kept.
   6. Everything else is removed.
4. **Clean** each genome with FCS-GX's own cleaning tool, applying only the removals.
5. **Combine** the cleaned genomes of a host into one reference and, for each index type listed in the host's `indexes` (none by default), build a bowtie2 or minimap2 index.

Every region, test result, rule and threshold is recorded, and the records can be written out as supplementary tables.

## Configuration

See [`config/config.example.yaml`](config/config.example.yaml).
A host is a name and its genomes (NCBI assembly accessions or local FASTA files).
Each genome's closest relative and the host's comparison genomes (reference ladder) are chosen from the NCBI taxonomy unless the config names them.
Each choice is saved in the library with every candidate considered, and reused on later runs while the genomes and settings are unchanged.
Every threshold is a setting with a default.

**Mitochondria.** Most nuclear assemblies are deposited without their mitochondrion; labs submit it as a separate NCBI Nucleotide record, or not at all. NCBI's RefSeq copies (GCF_) carry the species' RefSeq mitochondrial genome (an NC_ record) added after assembly. A genome entry can name that record, `{accession: GCA_059577405.1, mitochondrion: NC_064185.1}`, and it is added to the host reference right after the genome. Before anything is written, BoarderCollie checks that the record is a complete mitochondrial genome of the length NCBI states, that it is from the genome's own species, and that the assembly does not already hold a mitochondrion; any failure stops the run. The record is listed in `members.tsv` with the genome it belongs to (`member_of`), and its checks are saved beside it. It is not screened with FCS-GX: it is one curated 15-18 kb sequence of the host's own species.

## Running

```
mamba env create -f workflow/envs/boardercollie.yaml -p ./env
./env/bin/pip install --no-deps -e .
./env/bin/snakemake -s workflow/Snakefile --configfile config.yaml --cores 64 --resources mem_mb=740000
```

FCS-GX needs its database (about 520 GB) in memory, so screening runs one genome at a time on a large-memory node.
It is much faster when the database is on local storage: copy it once per run with `scripts/stage_fcs_gx_db.sh` and pass the copy with `--config fcs_gx_db=<copy>`, deleting the copy afterwards.

FCS-GX reports made earlier can be used instead of running FCS-GX.
Put them in one folder and name it in the config as `fcs_gx_reports`.
Each genome (and each relative or comparison genome) needs both its reports, named after its accession or FASTA file name: `<accession>.fcs_gx_report.txt` and `<accession>.taxonomy.rpt`, or as FCS-GX names them, `<accession>.<tax-id>.fcs_gx_report.txt` and `<accession>.<tax-id>.taxonomy.rpt`.
Before they are used, every sequence the reports name must be in the genome with the same length, and both reports must come from the same FCS-GX run; otherwise the run stops.
FCS-GX report headers do not state the tax-id it was run with, so that cannot be checked; the division FCS-GX assumed (for example `anml:primates`) is recorded in `check.json`.
Genomes without reports in the folder are screened as usual, which needs `fcs_gx_db`.
All reports must come from one database build: the build is read from the report headers, or given as `fcs_gx_db_build`.
Results made with these reports get their own run key, and `genome_summary.tsv` says where each genome's FCS-GX results came from.

Outputs, under `library_dir`:

- `genomes/<accession>/original/`: the genome as downloaded and its NCBI record (md5 checked against NCBI).
- `genomes/<accession>/fcs_gx/<database build>/`: FCS-GX reports.
- `genomes/<accession>/fcs_gx/precomputed/<reports hash>/`: copies of reports made earlier, and `check.json` (what was checked against the genome).
- `genomes/<accession>/runs/<run key>/`: `decisions.tsv`, the cleaned and removed sequence, the applied FCS-GX report and `run.json`.
- `hosts/<name>/<host key>/`: the joined host reference and `members.tsv`.

and under `index_dir`: `bt2/<name>/<host key>/` and `minimap2/<name>/<host key>/`.

A host with `reference_out: <path>` in the config also gets its joined reference at that fixed path (a hard link, or a copy on another file system), with `<path>.members.tsv` beside it.
Likewise `bowtie2_out: <stem>` links the host's bowtie2 index to `<stem>.1.bt2l` ... `<stem>.rev.2.bt2l` with `<stem>.index_done.txt`, and `minimap2_out: <path>` links its minimap2 index to `<path>` (each needs that index type in `indexes`).
The indexes themselves are built once per host key under `index_dir`; every workflow that asks for the same host links to them instead of building its own.

## Use from another Snakemake workflow

Another workflow (such as MAGmaker) can load BoarderCollie as a Snakemake module and ask for the cleaned host reference at a path it chooses:

```
module boardercollie:
    snakefile: "/path/to/BoarderCollie/workflow/Snakefile"
    config: config["boardercollie"]

use rule * from boardercollie as boardercollie_*
```

The `boardercollie` section of the calling workflow's config is a BoarderCollie config.
Give each host `reference_out` set to the path the calling workflow reads, for example MAGmaker's `host_filter: genome`; asking for that file runs everything needed to make it.
Set `bowtie2_out` and `minimap2_out` the same way to read the shared indexes at fixed paths; the first caller builds them, later callers only link.
Point `snakefile` at a local clone of this repository: the workflow takes the `boardercollie` package and `resources/` from that clone, so the package need not be installed in the calling workflow's environment.
Every step names `workflow/envs/boardercollie.yaml`, used when the calling workflow runs with `--use-conda`.
Snakemake records how each file was made in the `.snakemake` folder of the directory it runs from; a library made in one directory without `--use-conda` and later run from the same directory with it would be remade once, so keep one way of running per working directory.

## Development

```
python -m venv .venv
.venv/bin/pip install -e ".[test]"
.venv/bin/pytest
```

## License

MIT; see [LICENSE](LICENSE).

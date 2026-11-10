"""NCBI Datasets REST API (v2): genome downloads, assembly records, taxonomy."""

from __future__ import annotations

import hashlib
import io
import json
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"


def _get(url: str, tries: int = 5) -> bytes:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(10 * (i + 1))
    raise RuntimeError("unreachable")


def genome_report(accession: str) -> dict:
    d = json.loads(_get(f"{API}/genome/accession/{accession}/dataset_report"))
    reports = d.get("reports") or []
    if len(reports) != 1:
        raise ValueError(f"{accession}: expected one assembly record, got {len(reports)}")
    return reports[0]


def genomes_under(taxid: int, reference_only: bool = False) -> list[dict]:
    """Every current, non-atypical assembly record under an NCBI taxon (all pages)."""
    params = {"filters.assembly_version": "current", "filters.exclude_atypical": "true", "page_size": "1000"}
    if reference_only:
        params["filters.reference_only"] = "true"
    out, token = [], None
    while True:
        q = dict(params, **({"page_token": token} if token else {}))
        d = json.loads(_get(f"{API}/genome/taxon/{int(taxid)}/dataset_report?{urllib.parse.urlencode(q)}"))
        out.extend(d.get("reports") or [])
        token = d.get("next_page_token")
        if not token:
            return out


def download_genome(accession: str, out_fasta: Path) -> str:
    """Download the genome FASTA, check it against NCBI's md5, write it; returns the md5."""
    url = f"{API}/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA"
    z = zipfile.ZipFile(io.BytesIO(_get(url)))
    fna = [n for n in z.namelist() if n.endswith(".fna")]
    if len(fna) != 1:
        raise ValueError(f"{accession}: expected one .fna in the download, got {fna}")
    data = z.read(fna[0])
    md5 = hashlib.md5(data).hexdigest()
    sums = z.read("md5sum.txt").decode().split("\n") if "md5sum.txt" in z.namelist() else []
    want = [l.split()[0] for l in sums if l.strip().endswith(fna[0])]
    if not want or want[0] != md5:
        raise ValueError(f"{accession}: md5 {md5} does not match NCBI's md5sum.txt ({want})")
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_fasta.with_suffix(out_fasta.suffix + ".partial")
    tmp.write_bytes(data)
    tmp.rename(out_fasta)
    return md5


def ranks_of(taxids) -> dict[int, dict[str, int]]:
    """tax-id -> {rank: tax-id} from the NCBI classification."""
    return {t: v["ranks"] for t, v in taxonomy(taxids).items()}


def taxonomy(taxids) -> dict[int, dict]:
    """tax-id -> {"parents": set of ancestor tax-ids, "ranks": {rank: tax-id}, "name": str}."""
    out: dict[int, dict] = {}
    ids = sorted({int(t) for t in taxids})
    for i in range(0, len(ids), 100):
        chunk = ",".join(map(str, ids[i:i + 100]))
        d = json.loads(_get(f"{API}/taxonomy/taxon/{urllib.parse.quote(chunk)}/dataset_report"))
        for r in d.get("reports", []):
            t = r.get("taxonomy") or {}
            if "tax_id" not in t:
                continue
            ranks = {k: v.get("id") for k, v in (t.get("classification") or {}).items() if isinstance(v, dict)}
            out[int(t["tax_id"])] = {
                "parents": set(int(p) for p in t.get("parents") or []),
                "ranks": ranks,
                "name": (t.get("current_scientific_name") or {}).get("name", ""),
            }
    return out

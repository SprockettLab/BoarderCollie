"""Random access to large FASTA files without loading them into memory.

An index in samtools faidx format (name, length, offset, bases per line, bytes
per line) is built on first use and saved beside the work files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class FaiEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_bytes: int


def build_fai(fasta: Path, fai: Path) -> None:
    entries = []
    with open(fasta, "rb") as fh:
        name = None
        length = offset = line_bases = line_bytes = 0
        pos = 0
        last_short = False
        for raw in fh:
            if raw.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_bytes))
                name = raw[1:].split()[0].decode()
                length = line_bases = line_bytes = 0
                offset = pos + len(raw)
                last_short = False
            else:
                bases = len(raw.rstrip(b"\r\n"))
                if bases:
                    if last_short:
                        raise ValueError(f"{fasta}: {name} has lines of unequal length")
                    if line_bases == 0:
                        line_bases, line_bytes = bases, len(raw)
                    elif bases != line_bases:
                        last_short = bases < line_bases
                        if bases > line_bases:
                            raise ValueError(f"{fasta}: {name} has lines of unequal length")
                    length += bases
            pos += len(raw)
        if name is not None:
            entries.append((name, length, offset, line_bases, line_bytes))
    with open(fai, "w") as out:
        for e in entries:
            out.write("\t".join(map(str, e)) + "\n")


class IndexedFasta:
    def __init__(self, fasta: str | Path, fai: str | Path | None = None):
        self.fasta = Path(fasta)
        self.fai = Path(fai) if fai else Path(str(fasta) + ".fai")
        if not self.fai.exists():
            build_fai(self.fasta, self.fai)
        self.entries: dict[str, FaiEntry] = {}
        for line in open(self.fai):
            n, ln, off, lb, lby = line.rstrip("\n").split("\t")[:5]
            self.entries[n] = FaiEntry(n, int(ln), int(off), int(lb), int(lby))
        self._fh = open(self.fasta, "rb")

    def lengths(self) -> dict[str, int]:
        return {n: e.length for n, e in self.entries.items()}

    def fetch(self, name: str, start: int, end: int) -> str:
        """Bases start..end of sequence name, 1-based and inclusive, upper case."""
        e = self.entries[name]
        if not 1 <= start <= end <= e.length:
            raise ValueError(f"{name}:{start}-{end} outside 1-{e.length}")
        lb, lby = e.line_bases, e.line_bytes
        s0, e0 = start - 1, end
        first = e.offset + (s0 // lb) * lby + s0 % lb
        last = e.offset + ((e0 - 1) // lb) * lby + (e0 - 1) % lb
        self._fh.seek(first)
        raw = self._fh.read(last - first + 1)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()
        if len(seq) != end - start + 1:
            raise ValueError(f"{name}:{start}-{end}: read {len(seq)} bases")
        return seq

    def close(self) -> None:
        self._fh.close()


def read_fasta(path: str | Path) -> dict[str, str]:
    """Small FASTA files only: name -> sequence (upper case)."""
    seqs: dict[str, list[str]] = {}
    name = None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line and name is not None:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def write_fasta(records, path: str | Path, width: int = 80) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")

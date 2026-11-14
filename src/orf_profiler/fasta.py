"""Read transcript FASTA files and parse GENCODE header lines.

GENCODE transcript FASTA headers are pipe-separated. For protein-coding
transcripts they also carry the position of the coding sequence (CDS)
within the transcript, for example::

    >ENST00000641515.2|ENSG00000186092.7|...|OR4F5-201|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|

Coordinates in the header are 1-based and inclusive, and the CDS range
includes the stop codon. Long non-coding RNA transcripts have no CDS field.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, TextIO

from Bio.SeqIO.FastaIO import SimpleFastaParser


@dataclass(frozen=True)
class Transcript:
    """One transcript: identifiers, sequence, and CDS position if annotated."""

    transcript_id: str
    gene_id: str
    gene_name: str
    sequence: str
    cds_start: Optional[int] = None  # 1-based, inclusive
    cds_end: Optional[int] = None  # 1-based, inclusive (includes the stop codon)

    @property
    def has_cds(self) -> bool:
        return self.cds_start is not None and self.cds_end is not None


def open_text(path: str | Path) -> TextIO:
    """Open a plain or gzip-compressed text file for reading."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_header(header: str) -> tuple[str, str, str, Optional[int], Optional[int]]:
    """Parse a FASTA header into (transcript_id, gene_id, gene_name, cds_start, cds_end).

    GENCODE headers are split on ``|``. Any other header is accepted too: the
    first word becomes the transcript ID and the transcript is treated as
    having no annotated CDS.
    """
    header = header.strip()
    if "|" not in header:
        transcript_id = header.split()[0] if header else ""
        return transcript_id, "", "", None, None

    fields = header.split("|")
    transcript_id = fields[0]
    gene_id = fields[1] if len(fields) > 1 else ""
    gene_name = fields[5] if len(fields) > 5 else ""

    cds_start = cds_end = None
    for field in fields:
        if field.startswith("CDS:"):
            start_text, _, end_text = field[4:].partition("-")
            cds_start, cds_end = int(start_text), int(end_text)
            break
    return transcript_id, gene_id, gene_name, cds_start, cds_end


def clean_sequence(sequence: str) -> str:
    """Upper-case a sequence and convert RNA (U) to DNA (T) letters."""
    return sequence.upper().replace("U", "T")


def read_transcripts(path: str | Path) -> Iterator[Transcript]:
    """Yield every transcript in a FASTA file (plain or .gz)."""
    with open_text(path) as handle:
        for header, sequence in SimpleFastaParser(handle):
            transcript_id, gene_id, gene_name, cds_start, cds_end = parse_header(header)
            yield Transcript(
                transcript_id=transcript_id,
                gene_id=gene_id,
                gene_name=gene_name,
                sequence=clean_sequence(sequence),
                cds_start=cds_start,
                cds_end=cds_end,
            )

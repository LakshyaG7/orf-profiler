"""Run the full analysis over a FASTA file: find, classify and characterize ORFs."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from itertools import islice
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Iterator, NamedTuple, Optional, TextIO

from . import __version__
from .classify import ANNOTATED, EXTENSION, classify
from .fasta import Transcript, read_transcripts
from .features import kozak_context, kozak_strength, peptide_features, translate
from .orfs import STOP_CODONS, find_orfs

COLUMNS = [
    "orf_id",
    "transcript_id",
    "gene_id",
    "gene_name",
    "orf_class",
    "frame",
    "start",
    "end",
    "start_codon",
    "stop_codon",
    "kozak_context",
    "kozak_strength",
    "aa_length",
    "molecular_weight",
    "isoelectric_point",
    "gravy",
]

BATCH_SIZE = 200  # transcripts handed to a worker process at a time


class OrfRecord(NamedTuple):
    """One output row. ``start`` and ``end`` are 1-based and inclusive."""

    orf_id: str
    transcript_id: str
    gene_id: str
    gene_name: str
    orf_class: str
    frame: int
    start: int
    end: int
    start_codon: str
    stop_codon: str
    kozak_context: str
    kozak_strength: str
    aa_length: int
    molecular_weight: Optional[float]
    isoelectric_point: Optional[float]
    gravy: Optional[float]
    peptide: str

    def to_row(self) -> str:
        def number(value: Optional[float], digits: int) -> str:
            return "NA" if value is None else f"{value:.{digits}f}"

        return "\t".join(
            [
                self.orf_id,
                self.transcript_id,
                self.gene_id,
                self.gene_name,
                self.orf_class,
                str(self.frame),
                str(self.start),
                str(self.end),
                self.start_codon,
                self.stop_codon,
                self.kozak_context,
                self.kozak_strength,
                str(self.aa_length),
                number(self.molecular_weight, 1),
                number(self.isoelectric_point, 2),
                number(self.gravy, 3),
            ]
        )


def profile_transcript(
    transcript: Transcript,
    min_aa: int = 10,
    max_aa: Optional[int] = None,
    near_cognate: bool = False,
) -> list[OrfRecord]:
    """Find every ORF in one transcript and describe it."""
    sequence = transcript.sequence
    records = []
    for orf in find_orfs(sequence, min_aa=min_aa, near_cognate=near_cognate):
        if max_aa is not None and orf.aa_length > max_aa:
            continue
        start, end = orf.start + 1, orf.end  # to 1-based, inclusive
        peptide = translate(sequence[orf.start : orf.end])
        features = peptide_features(peptide)
        records.append(
            OrfRecord(
                orf_id=f"{transcript.transcript_id}:{start}-{end}",
                transcript_id=transcript.transcript_id,
                gene_id=transcript.gene_id,
                gene_name=transcript.gene_name,
                orf_class=classify(start, end, transcript.cds_start, transcript.cds_end),
                frame=orf.frame + 1,
                start=start,
                end=end,
                start_codon=orf.start_codon,
                stop_codon=orf.stop_codon,
                kozak_context=kozak_context(sequence, orf.start),
                kozak_strength=kozak_strength(sequence, orf.start),
                aa_length=orf.aa_length,
                molecular_weight=features.molecular_weight if features else None,
                isoelectric_point=features.isoelectric_point if features else None,
                gravy=features.gravy if features else None,
                peptide=peptide,
            )
        )
    return records


def cds_is_complete(transcript: Transcript) -> bool:
    """True if the annotated CDS is itself a complete ATG-to-stop ORF.

    GENCODE includes transcripts whose CDS is incomplete at either end, starts
    at a non-ATG codon, or reads through a stop codon (selenoproteins). Those
    cannot be recovered by an ORF scan and are excluded from the sanity check.
    """
    if not transcript.has_cds:
        return False
    cds = transcript.sequence[transcript.cds_start - 1 : transcript.cds_end]
    if len(cds) < 6 or len(cds) % 3 != 0:
        return False
    if not cds.startswith("ATG") or cds[-3:] not in STOP_CODONS:
        return False
    return all(cds[i : i + 3] not in STOP_CODONS for i in range(0, len(cds) - 3, 3))


def _process_batch(args) -> tuple[str, str, Counter]:
    """Worker: profile a batch of transcripts and return text chunks plus counts."""
    transcripts, min_aa, max_aa, near_cognate = args
    rows, fasta, stats = [], [], Counter()
    for transcript in transcripts:
        records = profile_transcript(transcript, min_aa, max_aa, near_cognate)
        stats["transcripts"] += 1
        stats["orfs"] += len(records)
        classes = set()
        for record in records:
            rows.append(record.to_row())
            fasta.append(f">{record.orf_id} {record.orf_class} {record.gene_name}\n{record.peptide}")
            classes.add(record.orf_class)
            stats["class", record.orf_class] += 1
            stats["length", record.orf_class, record.aa_length] += 1
            stats["kozak", record.orf_class, record.kozak_strength] += 1
            stats["start_codon", record.orf_class, record.start_codon] += 1

        if transcript.has_cds:
            stats["transcripts_with_cds"] += 1
            # The check below only makes sense if the CDS itself passes the length filter.
            cds_aa = (transcript.cds_end - transcript.cds_start + 1) // 3 - 1
            in_range = cds_aa >= min_aa and (max_aa is None or cds_aa <= max_aa)
            if in_range and cds_is_complete(transcript):
                stats["cds_complete"] += 1
                if ANNOTATED in classes:
                    stats["cds_recovered_exact"] += 1
                elif EXTENSION in classes:
                    stats["cds_recovered_as_extension"] += 1
    return "\n".join(rows), "\n".join(fasta), stats


def _batches(items: Iterable, size: int) -> Iterator[list]:
    iterator = iter(items)
    while batch := list(islice(iterator, size)):
        yield batch


def _open_out(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "wt", compresslevel=4)
    return open(path, "w")


def scan_fasta(
    fasta_path: str | Path,
    out_dir: str | Path,
    min_aa: int = 10,
    max_aa: Optional[int] = None,
    near_cognate: bool = False,
    threads: int = 1,
    compress: bool = False,
    write_peptides: bool = True,
) -> dict:
    """Scan every transcript in ``fasta_path`` and write results to ``out_dir``.

    Writes ``orfs.tsv`` (one row per ORF), ``peptides.fa`` and the summary
    tables used for plotting. Returns the run summary as a dictionary.
    """
    from .report import write_summaries  # imported here to keep worker start-up light

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".gz" if compress else ""
    tsv_path = out_dir / f"orfs.tsv{suffix}"
    fasta_out_path = out_dir / f"peptides.fa{suffix}"

    jobs = (
        (batch, min_aa, max_aa, near_cognate)
        for batch in _batches(read_transcripts(fasta_path), BATCH_SIZE)
    )
    stats: Counter = Counter()

    tsv = _open_out(tsv_path)
    peptides = _open_out(fasta_out_path) if write_peptides else None
    pool = Pool(threads) if threads > 1 else None
    try:
        tsv.write("\t".join(COLUMNS) + "\n")
        # Hand the pool a few batches at a time so the whole file is never held in memory.
        for group in _batches(jobs, threads * 4):
            results = pool.map(_process_batch, group) if pool else map(_process_batch, group)
            for rows, fasta, batch_stats in results:
                if rows:
                    tsv.write(rows + "\n")
                    if peptides:
                        peptides.write(fasta + "\n")
                stats.update(batch_stats)
    finally:
        if pool:
            pool.close()
            pool.join()
        tsv.close()
        if peptides:
            peptides.close()

    run_info = {
        "orf_profiler_version": __version__,
        "input": str(fasta_path),
        "parameters": {"min_aa": min_aa, "max_aa": max_aa, "near_cognate": near_cognate},
        "transcripts": stats["transcripts"],
        "transcripts_with_cds": stats["transcripts_with_cds"],
        "orfs": stats["orfs"],
        "cds_check": {
            "complete_cds": stats["cds_complete"],
            "recovered_exact": stats["cds_recovered_exact"],
            "recovered_as_extension": stats["cds_recovered_as_extension"],
        },
    }
    write_summaries(out_dir, stats)
    with open(out_dir / "run_info.json", "w") as handle:
        json.dump(run_info, handle, indent=2)
        handle.write("\n")
    return run_info

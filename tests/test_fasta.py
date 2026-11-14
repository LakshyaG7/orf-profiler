import gzip
from pathlib import Path

from orf_profiler.fasta import clean_sequence, parse_header, read_transcripts

DATA = Path(__file__).parent / "data"

PC_HEADER = (
    "ENST00000641515.2|ENSG00000186092.7|OTTHUMG00000001094.4|OTTHUMT00000003223.4|"
    "OR4F5-201|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|"
)
LNC_HEADER = "ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|"


def test_parse_protein_coding_header():
    assert parse_header(PC_HEADER) == ("ENST00000641515.2", "ENSG00000186092.7", "OR4F5", 61, 1041)


def test_parse_lncrna_header_has_no_cds():
    assert parse_header(LNC_HEADER) == ("ENST00000832824.1", "ENSG00000290825.2", "DDX11L16", None, None)


def test_parse_header_without_utrs():
    # A transcript with no annotated UTRs still has a CDS field.
    header = "ENST1.1|ENSG1.1|-|-|GENE-201|GENE|300|CDS:1-300|"
    assert parse_header(header)[3:] == (1, 300)


def test_parse_plain_header():
    assert parse_header("my_transcript some description") == ("my_transcript", "", "", None, None)


def test_clean_sequence_handles_lowercase_and_rna():
    assert clean_sequence("augcUa") == "ATGCTA"


def test_read_transcripts():
    transcripts = list(read_transcripts(DATA / "tiny.fa"))
    assert len(transcripts) == 5
    smim1 = transcripts[0]
    assert smim1.gene_name == "SMIM1"
    assert len(smim1.sequence) == 506
    assert (smim1.cds_start, smim1.cds_end) == (235, 471)
    assert smim1.has_cds
    assert not transcripts[-1].has_cds


def test_read_gzipped_fasta(tmp_path):
    path = tmp_path / "tiny.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write((DATA / "tiny.fa").read_text())
    assert [t.transcript_id for t in read_transcripts(path)] == [
        t.transcript_id for t in read_transcripts(DATA / "tiny.fa")
    ]

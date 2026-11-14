import csv
import gzip
import json
from pathlib import Path

from orf_profiler.cli import main
from orf_profiler.fasta import Transcript, read_transcripts
from orf_profiler.scan import COLUMNS, cds_is_complete, profile_transcript, scan_fasta

DATA = Path(__file__).parent / "data"
SMIM1_PROTEIN = "MQPQESHVHYSRWEDGSRDGVSLGAVSSTEEASRCRRISQRLCTGKLGIAMKVLGGVALFWIIFILGYLTGYYVHKCK"


def read_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def smim1():
    return next(read_transcripts(DATA / "tiny.fa"))


def test_annotated_cds_of_a_real_microprotein_is_recovered():
    records = profile_transcript(smim1())
    annotated = [r for r in records if r.orf_class == "annotated"]
    assert len(annotated) == 1
    record = annotated[0]
    assert (record.start, record.end) == (235, 471)
    assert record.aa_length == 78
    assert record.peptide == SMIM1_PROTEIN
    assert record.orf_id == "ENST00000642557.4:235-471"


def test_coordinates_are_one_based_inclusive():
    transcript = smim1()
    for record in profile_transcript(transcript, near_cognate=True):
        orf_sequence = transcript.sequence[record.start - 1 : record.end]
        assert orf_sequence[:3] == record.start_codon
        assert orf_sequence[-3:] == record.stop_codon
        assert len(orf_sequence) == 3 * (record.aa_length + 1)


def test_near_cognate_adds_orfs_without_removing_atg_orfs():
    transcript = smim1()
    atg_only = {r.orf_id for r in profile_transcript(transcript)}
    with_near = {r.orf_id for r in profile_transcript(transcript, near_cognate=True)}
    assert atg_only < with_near


def test_max_aa_filter():
    records = profile_transcript(smim1(), near_cognate=True, max_aa=50)
    assert records and all(r.aa_length <= 50 for r in records)


def test_lncrna_transcript_orfs():
    lnc = list(read_transcripts(DATA / "tiny.fa"))[-1]
    records = profile_transcript(lnc)
    assert [r.orf_class for r in records] == ["lncRNA-ORF"]
    assert records[0].kozak_strength == "strong"  # GCCACC ATG G


def test_cds_is_complete():
    good = Transcript("t", "g", "n", "CC" + "ATGGCTGCTTAA" + "CC", 3, 14)
    assert cds_is_complete(good)
    no_stop = Transcript("t", "g", "n", "CC" + "ATGGCTGCTGCT" + "CC", 3, 14)
    assert not cds_is_complete(no_stop)
    non_atg = Transcript("t", "g", "n", "CC" + "CTGGCTGCTTAA" + "CC", 3, 14)
    assert not cds_is_complete(non_atg)
    internal_stop = Transcript("t", "g", "n", "CC" + "ATGTGAGCTTAA" + "CC", 3, 14)
    assert not cds_is_complete(internal_stop)
    assert not cds_is_complete(Transcript("t", "g", "n", "ATGGCTTAA"))


def test_scan_writes_all_outputs(tmp_path):
    info = scan_fasta(DATA / "tiny.fa", tmp_path, near_cognate=True)
    rows = read_rows(tmp_path / "orfs.tsv")
    assert list(rows[0]) == COLUMNS
    assert len(rows) == info["orfs"] > 0
    assert info["transcripts"] == 5
    assert info["transcripts_with_cds"] == 4
    assert info["cds_check"] == {"complete_cds": 4, "recovered_exact": 4, "recovered_as_extension": 0}

    peptides = (tmp_path / "peptides.fa").read_text()
    assert peptides.count(">") == len(rows)
    assert SMIM1_PROTEIN in peptides

    summary = {r["orf_class"]: int(r["n_orfs"]) for r in read_rows(tmp_path / "class_summary.tsv")}
    assert sum(summary.values()) == len(rows)
    assert summary["annotated"] == 4
    assert json.loads((tmp_path / "run_info.json").read_text())["parameters"]["near_cognate"] is True


def test_scan_gives_same_result_with_several_processes(tmp_path):
    scan_fasta(DATA / "tiny.fa", tmp_path / "one", near_cognate=True, threads=1)
    scan_fasta(DATA / "tiny.fa", tmp_path / "two", near_cognate=True, threads=2)
    assert (tmp_path / "one" / "orfs.tsv").read_text() == (tmp_path / "two" / "orfs.tsv").read_text()


def test_command_line_scan_and_plot(tmp_path, capsys):
    results = tmp_path / "results"
    assert main(["scan", str(DATA / "tiny.fa"), "-o", str(results), "--near-cognate", "--gzip"]) == 0
    assert "Found" in capsys.readouterr().out
    assert read_rows(results / "orfs.tsv.gz")

    figures = tmp_path / "figures"
    assert main(["plot", str(results), "-o", str(figures)]) == 0
    assert (figures / "orf_length_by_class.png").stat().st_size > 0
    assert (figures / "kozak_by_class.png").stat().st_size > 0


def test_command_line_rejects_bad_length_limits(tmp_path):
    assert main(["scan", str(DATA / "tiny.fa"), "-o", str(tmp_path), "--min-aa", "50", "--max-aa", "20"]) == 2

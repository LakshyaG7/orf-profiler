import pytest

from orf_profiler.classify import classify

# CDS occupies positions 101-400 (1-based, inclusive, including its stop codon).
CDS = (101, 400)


@pytest.mark.parametrize(
    "orf, expected",
    [
        ((101, 400), "annotated"),
        ((41, 400), "extension"),  # same stop, starts upstream
        ((161, 400), "truncation"),  # same stop, starts downstream
        ((10, 60), "uORF"),  # entirely in the 5' UTR
        ((10, 100), "uORF"),  # ends on the last base before the CDS
        ((50, 150), "uoORF"),  # starts in the 5' UTR, ends inside the CDS
        ((50, 450), "uoORF"),  # starts upstream and runs past the CDS
        ((150, 250), "intORF"),  # inside the CDS, different frame
        ((350, 450), "doORF"),  # starts in the CDS, ends in the 3' UTR
        ((401, 460), "dORF"),  # starts on the first base after the CDS
        ((420, 480), "dORF"),
    ],
)
def test_classes_relative_to_cds(orf, expected):
    assert classify(*orf, *CDS) == expected


def test_no_cds_is_lncrna_orf():
    assert classify(10, 60, None, None) == "lncRNA-ORF"

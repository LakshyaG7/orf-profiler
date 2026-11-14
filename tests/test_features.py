import pytest

from orf_profiler.features import kozak_context, kozak_strength, peptide_features, translate


def test_translate_drops_stop_codon():
    assert translate("ATGGCTAAATAA") == "MAK"


def test_translate_near_cognate_start_is_methionine():
    # CTG codes for leucine internally but initiates with methionine.
    assert translate("CTGCTGTAA") == "ML"


def test_translate_ambiguous_codon():
    assert translate("ATGGNTTAA") == "MX"


#           -6 ..... -1  +1+2+3 +4
STRONG = "GCCACC" + "ATG" + "G"  # A at -3, G at +4
ONLY_MINUS3 = "GCCGCC" + "ATG" + "C"
ONLY_PLUS4 = "GCCTCC" + "ATG" + "G"
WEAK = "GCCTCC" + "ATG" + "C"


@pytest.mark.parametrize(
    "sequence, expected",
    [(STRONG, "strong"), (ONLY_MINUS3, "moderate"), (ONLY_PLUS4, "moderate"), (WEAK, "weak")],
)
def test_kozak_strength(sequence, expected):
    assert kozak_strength(sequence, 6) == expected


def test_kozak_unknown_at_transcript_ends():
    assert kozak_strength("ATGGCC", 0) == "NA"  # no -3 position
    assert kozak_strength("GCCACCATG", 6) == "NA"  # no +4 position


def test_kozak_context_window():
    assert kozak_context(STRONG + "TTT", 6) == "GCCACCATGG"
    assert kozak_context("ATGGCC", 0) == "ATGG"


def test_peptide_features_match_known_values():
    # Human insulin B chain (30 aa); reference values from Biopython's ProtParam.
    features = peptide_features("FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    assert features.molecular_weight == pytest.approx(3429.9, abs=0.5)
    assert features.isoelectric_point == pytest.approx(6.90, abs=0.05)
    assert features.gravy == pytest.approx(0.220, abs=0.005)


def test_hydrophobic_peptide_has_positive_gravy():
    assert peptide_features("MLLVVIIAALL").gravy > 0
    assert peptide_features("MDDEEKKRRDD").gravy < 0


def test_peptide_with_unknown_residue_has_no_features():
    assert peptide_features("MAXK") is None
    assert peptide_features("") is None

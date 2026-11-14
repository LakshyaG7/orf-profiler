from orf_profiler.orfs import find_orfs

# A 10-codon ORF: ATG + 9 codons + stop = 10 amino acids.
BODY = "GCT" * 9
ORF10 = "ATG" + BODY + "TAA"


def orfs(sequence, **kwargs):
    return list(find_orfs(sequence, **kwargs))


def test_single_orf_in_frame_0():
    found = orfs(ORF10)
    assert len(found) == 1
    orf = found[0]
    assert (orf.start, orf.end, orf.frame) == (0, 33, 0)
    assert (orf.start_codon, orf.stop_codon) == ("ATG", "TAA")
    assert orf.aa_length == 10


def test_orf_found_in_each_frame():
    for shift in (1, 2):
        found = orfs("C" * shift + ORF10)
        assert len(found) == 1
        assert found[0].frame == shift
        assert found[0].start == shift


def test_all_three_stop_codons():
    for stop in ("TAA", "TAG", "TGA"):
        found = orfs("ATG" + BODY + stop)
        assert [o.stop_codon for o in found] == [stop]


def test_no_stop_codon_means_no_orf():
    assert orfs("ATG" + BODY + "GCT") == []


def test_no_start_codon_means_no_orf():
    assert orfs("GCT" * 12 + "TAA") == []


def test_longest_orf_per_stop_codon():
    # Two in-frame ATGs before one stop: only the upstream one is reported.
    sequence = "ATG" + "GCT" * 5 + "ATG" + "GCT" * 9 + "TAA"
    found = orfs(sequence)
    assert len(found) == 1
    assert found[0].start == 0


def test_stop_codon_resets_the_search():
    # Short ORF, stop, then a second ORF in the same frame.
    sequence = "ATG" + "GCT" * 2 + "TAA" + ORF10
    found = orfs(sequence)
    assert len(found) == 1
    assert found[0].start == 12


def test_min_aa_filter():
    assert orfs(ORF10, min_aa=11) == []
    assert len(orfs(ORF10, min_aa=10)) == 1


def test_out_of_frame_stop_is_ignored():
    # "TAA" appears at index 4 but in frame 1, so the frame-0 ORF is intact.
    sequence = "ATGCTAAGC" + "GCT" * 8 + "TAG"
    found = orfs(sequence)
    assert [o.frame for o in found] == [0]
    assert found[0].end == len(sequence)


def test_near_cognate_off_by_default():
    assert orfs("CTG" + BODY + "TAA") == []


def test_near_cognate_start_without_atg():
    for start in ("CTG", "GTG", "TTG"):
        found = orfs(start + BODY + "TAA", near_cognate=True)
        assert len(found) == 1
        assert found[0].start_codon == start
        assert found[0].aa_length == 10


def test_near_cognate_upstream_of_atg_gives_two_orfs():
    sequence = "CTG" + "GCT" * 3 + ORF10
    found = orfs(sequence, near_cognate=True)
    assert [(o.start, o.start_codon) for o in found] == [(0, "CTG"), (12, "ATG")]
    assert found[0].end == found[1].end


def test_near_cognate_downstream_of_atg_is_not_reported():
    sequence = "ATG" + "GCT" * 3 + "CTG" + "GCT" * 9 + "TAA"
    found = orfs(sequence, near_cognate=True)
    assert [(o.start, o.start_codon) for o in found] == [(0, "ATG")]


def test_sequence_shorter_than_a_codon():
    assert orfs("AT") == []
    assert orfs("") == []

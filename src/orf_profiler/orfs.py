"""Find open reading frames (ORFs) in a transcript sequence.

A transcript is read in the 5' to 3' direction, so only the three forward
reading frames are scanned. An ORF runs from a start codon to the first
in-frame stop codon. ORFs with no stop codon before the end of the
transcript are not reported.
"""

from __future__ import annotations

from typing import Iterator, NamedTuple

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
CANONICAL_START = "ATG"
# Near-cognate start codons differ from ATG by one base. CTG, GTG and TTG are
# common alternative initiation codons; others (such as ACG) are not searched.
NEAR_COGNATE_STARTS = frozenset({"CTG", "GTG", "TTG"})


class ORF(NamedTuple):
    """An ORF in 0-based, half-open coordinates on the transcript.

    ``start`` is the index of the first base of the start codon and ``end`` is
    one past the last base of the stop codon, so ``sequence[start:end]`` is the
    whole ORF including its stop codon.
    """

    start: int
    end: int
    frame: int
    start_codon: str
    stop_codon: str

    @property
    def aa_length(self) -> int:
        """Number of amino acids encoded (the stop codon is not counted)."""
        return (self.end - self.start) // 3 - 1


def find_orfs(sequence: str, min_aa: int = 10, near_cognate: bool = False) -> Iterator[ORF]:
    """Yield the ORFs in ``sequence`` (upper-case DNA letters).

    For every stop codon, the longest ORF ending there is reported, i.e. the
    one that begins at the most upstream ATG since the previous in-frame stop.

    With ``near_cognate=True``, a second ORF is reported for the same stop
    codon when a CTG, GTG or TTG lies upstream of the first ATG (or when there
    is no ATG at all). It represents a possible non-canonical start.

    ORFs encoding fewer than ``min_aa`` amino acids are skipped.
    """
    length = len(sequence)
    for frame in range(3):
        first_atg = None  # position of the first ATG since the last stop
        first_near = None  # position of the first near-cognate start before any ATG
        for i in range(frame, length - 2, 3):
            codon = sequence[i : i + 3]
            if codon in STOP_CODONS:
                end = i + 3
                if first_near is not None and (i - first_near) // 3 >= min_aa:
                    yield ORF(first_near, end, frame, sequence[first_near : first_near + 3], codon)
                if first_atg is not None and (i - first_atg) // 3 >= min_aa:
                    yield ORF(first_atg, end, frame, CANONICAL_START, codon)
                first_atg = None
                first_near = None
            elif first_atg is None:
                if codon == CANONICAL_START:
                    first_atg = i
                elif near_cognate and first_near is None and codon in NEAR_COGNATE_STARTS:
                    first_near = i

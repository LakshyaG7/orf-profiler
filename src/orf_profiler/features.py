"""Sequence features of an ORF and of the peptide it encodes."""

from __future__ import annotations

from functools import lru_cache
from typing import NamedTuple, Optional

from Bio.Data import CodonTable
from Bio.SeqUtils.ProtParam import ProteinAnalysis

_CODON_TABLE = CodonTable.unambiguous_dna_by_id[1].forward_table  # standard genetic code
_STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

KOZAK_STRONG = "strong"
KOZAK_MODERATE = "moderate"
KOZAK_WEAK = "weak"
KOZAK_UNKNOWN = "NA"
KOZAK_ORDER = [KOZAK_STRONG, KOZAK_MODERATE, KOZAK_WEAK, KOZAK_UNKNOWN]


def translate(orf_sequence: str) -> str:
    """Translate an ORF (start codon through stop codon) into a peptide.

    The stop codon is dropped. The first residue is always methionine:
    initiation uses the initiator Met-tRNA even at near-cognate start codons.
    Codons containing an ambiguous base (N) become ``X``.
    """
    residues = ["M"]
    for i in range(3, len(orf_sequence) - 3, 3):
        residues.append(_CODON_TABLE.get(orf_sequence[i : i + 3], "X"))
    return "".join(residues)


def kozak_context(sequence: str, start: int) -> str:
    """Return the bases from -6 to +4 around a start codon (0-based ``start``).

    Positions are numbered with the first base of the start codon as +1, so
    the returned string covers ``-6 .. -1``, the start codon, and ``+4``. It is
    shorter if the start codon is too close to either end of the transcript.
    """
    return sequence[max(0, start - 6) : start + 4]


def kozak_strength(sequence: str, start: int) -> str:
    """Score the Kozak context of a start codon as strong, moderate or weak.

    The two positions that matter most for start-codon recognition are -3
    (a purine, A or G) and +4 (a G). Both present is ``strong``, one is
    ``moderate`` and neither is ``weak``. If either position falls outside the
    transcript the result is ``NA``.
    """
    minus3 = start - 3
    plus4 = start + 3
    if minus3 < 0 or plus4 >= len(sequence):
        return KOZAK_UNKNOWN
    matches = (sequence[minus3] in "AG") + (sequence[plus4] == "G")
    return (KOZAK_WEAK, KOZAK_MODERATE, KOZAK_STRONG)[matches]


class PeptideFeatures(NamedTuple):
    molecular_weight: float  # Daltons
    isoelectric_point: float  # pH at which the net charge is zero
    gravy: float  # grand average of hydropathy (Kyte-Doolittle); > 0 is hydrophobic


@lru_cache(maxsize=200_000)
def peptide_features(peptide: str) -> Optional[PeptideFeatures]:
    """Compute physicochemical features of a peptide with Biopython.

    Returns None for peptides containing a non-standard residue (for example
    ``X`` from an ambiguous codon), because the features are undefined.

    Results are cached: transcript isoforms of one gene share most of their
    ORFs, so the same peptide is seen many times in a GENCODE-scale run.
    """
    if not peptide or not _STANDARD_AMINO_ACIDS.issuperset(peptide):
        return None
    analysis = ProteinAnalysis(peptide)
    return PeptideFeatures(
        molecular_weight=analysis.molecular_weight(),
        isoelectric_point=analysis.isoelectric_point(),
        gravy=analysis.gravy(),
    )

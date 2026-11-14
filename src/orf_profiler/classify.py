"""Classify an ORF by where it sits relative to the annotated coding sequence.

All coordinates here are 1-based and inclusive, matching GENCODE headers, and
both the ORF end and the CDS end include the stop codon.

Classes (names follow the convention used for ribosome-profiling ORFs):

==============  ==============================================================
annotated       same start and stop as the annotated CDS
extension       same stop as the CDS, starts upstream (N-terminal extension)
truncation      same stop as the CDS, starts downstream of the annotated start
uORF            upstream ORF: ends before the CDS starts (inside the 5' UTR)
uoORF           upstream overlapping ORF: starts in the 5' UTR, runs into the CDS
intORF          internal ORF: lies inside the CDS in a different reading frame
doORF           downstream overlapping ORF: starts in the CDS, ends in the 3' UTR
dORF            downstream ORF: starts after the CDS ends (inside the 3' UTR)
lncRNA-ORF      the transcript has no annotated CDS
==============  ==============================================================
"""

from __future__ import annotations

from typing import Optional

ANNOTATED = "annotated"
EXTENSION = "extension"
TRUNCATION = "truncation"
UORF = "uORF"
UOORF = "uoORF"
INTORF = "intORF"
DOORF = "doORF"
DORF = "dORF"
LNCRNA_ORF = "lncRNA-ORF"

# Display order, 5' to 3' along the transcript, with lncRNA ORFs last.
CLASS_ORDER = [ANNOTATED, EXTENSION, TRUNCATION, UORF, UOORF, INTORF, DOORF, DORF, LNCRNA_ORF]


def classify(
    orf_start: int,
    orf_end: int,
    cds_start: Optional[int],
    cds_end: Optional[int],
) -> str:
    """Return the class of an ORF given the CDS position (or None if no CDS)."""
    if cds_start is None or cds_end is None:
        return LNCRNA_ORF

    # Two ORFs that end at the same stop codon are in the same reading frame.
    if orf_end == cds_end:
        if orf_start == cds_start:
            return ANNOTATED
        return EXTENSION if orf_start < cds_start else TRUNCATION

    if orf_end < cds_start:
        return UORF
    if orf_start > cds_end:
        return DORF
    if orf_start < cds_start:
        return UOORF
    if orf_end > cds_end:
        return DOORF
    return INTORF

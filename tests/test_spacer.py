"""Tests for spacer extraction."""

from scalpel.design import SpacerExtractor, count_pam_sites
from scalpel.models.enums import CasVariantType, Strand


class TestSpacerExtraction:
    """Tests for spacer extraction from sequences."""

    def test_extract_spacers_both_strands(self) -> None:
        seq = "ATCGATCGATCGATCGATCGAGGATCGATCGATCGATCGATCGAGG"
        spacers = SpacerExtractor(CasVariantType.SPCAS9).extract_spacers(seq)
        assert spacers
        strands = {s.strand for s in spacers}
        assert Strand.PLUS in strands
        assert Strand.MINUS in strands

    def test_spacer_length_matches_spcas9(self) -> None:
        seq = "ATCGATCGATCGATCGATCGAGGATCGATCGATCGATCGATCGAGG"
        spacers = SpacerExtractor(CasVariantType.SPCAS9).extract_spacers(seq)
        assert spacers
        assert all(len(s.spacer_sequence) == 20 for s in spacers)

    def test_spacer_length_matches_cas12a(self) -> None:
        seq = "TTTA" + ("A" * 23) + "GCGCGC" + "TTTG" + ("C" * 23)
        spacers = SpacerExtractor(CasVariantType.CAS12A).extract_spacers(seq)
        assert spacers
        assert all(len(s.spacer_sequence) == 23 for s in spacers)

    def test_cut_site_calculation_cas9(self) -> None:
        # 20bp spacer + AGG PAM; cut is 3bp upstream of PAM start.
        seq = "ATCGATCGATCGATCGATCGAGG"
        spacers = SpacerExtractor(CasVariantType.SPCAS9).extract_spacers(
            seq,
            chromosome="chr1",
            start_position=100,
        )
        assert spacers
        plus = [s for s in spacers if s.strand == Strand.PLUS]
        assert plus
        assert plus[0].cut_site == 117

    def test_no_spacers_without_pam(self) -> None:
        no_pam_seq = "ATATATATATATATATATATATATATATATATATATATATATATATATATAT"
        spacers = SpacerExtractor(CasVariantType.SPCAS9).extract_spacers(no_pam_seq)
        assert spacers == []

    def test_invalid_spacers_with_n_are_filtered(self) -> None:
        # Only candidate spacer contains N, so candidate should be rejected.
        seq = ("N" * 20) + "AGG"
        spacers = SpacerExtractor(CasVariantType.SPCAS9).extract_spacers(seq)
        assert spacers == []

    def test_count_pam_sites_consistency(self) -> None:
        seq = "ATCGATCGATCGATCGATCGAGGATCGATCGATCGATCGATCGAGG"
        estimated = count_pam_sites(seq, CasVariantType.SPCAS9)
        extracted = len(SpacerExtractor(CasVariantType.SPCAS9).extract_spacers(seq))
        # Extraction can include reverse-strand matches beyond forward PAM counting.
        assert extracted >= estimated

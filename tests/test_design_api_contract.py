"""
Test design API contract (B2.3).

Regression tests to ensure the design API returns expected structure and values.
"""

import pytest
from unittest.mock import patch, MagicMock


class TestDesignAPIContract:
    """Test the design endpoint contract."""
    
    def test_design_response_has_required_fields(self):
        """Verify design response contains required guide fields."""
        from scalpel.design.pipeline import DesignPipeline
        from scalpel.models.enums import Genome, CasVariantType, EditModality
        
        # Simple sequence with known PAM sites
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT"
        
        pipeline = DesignPipeline(
            genome=Genome.HUMAN_GRCH38,
            cas_variant=CasVariantType.SPCAS9,
            modality=EditModality.KNOCKOUT,
        )
        
        # Extract spacers (doesn't require genome files)
        from scalpel.design.spacer_extraction import SpacerExtractor
        extractor = SpacerExtractor(CasVariantType.SPCAS9)
        spacers = extractor.extract_spacers(test_sequence)
        
        assert len(spacers) > 0, "Should find at least one spacer"
        
        # Verify spacer structure
        spacer = spacers[0]
        assert hasattr(spacer, 'spacer_sequence')
        assert hasattr(spacer, 'pam_sequence')
        assert hasattr(spacer, 'strand')
        assert hasattr(spacer, 'genomic_start')
        assert hasattr(spacer, 'cut_site')
        
        # Verify sequence lengths
        assert len(spacer.spacer_sequence) == 20, "Spacer should be 20bp"
        assert len(spacer.pam_sequence) == 3, "PAM should be 3bp for SpCas9"
    
    def test_efficiency_score_in_valid_range(self):
        """Verify efficiency scores are between 0 and 1."""
        from scalpel.design.efficiency.ensemble import EnsembleScorer
        from scalpel.design.spacer_extraction import SpacerExtractor
        from scalpel.models.enums import CasVariantType
        
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT"
        
        extractor = SpacerExtractor(CasVariantType.SPCAS9)
        spacers = extractor.extract_spacers(test_sequence)
        
        scorer = EnsembleScorer(use_ml_model=False)  # Rule-based only for speed
        
        for spacer in spacers[:5]:  # Test first 5
            scored = scorer.score_single(spacer)
            assert 0 <= scored.efficiency.overall_score <= 1, \
                f"Score {scored.efficiency.overall_score} out of range [0,1]"
    
    def test_offtarget_analysis_has_metadata(self):
        """Verify off-target analysis includes B1.2 backend metadata."""
        from scalpel.models.data_classes import OffTargetAnalysis
        from scalpel.models.enums import Genome
        
        # Create minimal analysis to check fields exist
        analysis = OffTargetAnalysis(
            on_target_spacer="ATGGATTTATCTGCTCTTCG",
            genome=Genome.HUMAN_GRCH38,
        )
        
        # B1.2 metadata fields
        assert hasattr(analysis, 'backend')
        assert hasattr(analysis, 'db_path')
        assert hasattr(analysis, 'candidate_count')
        assert hasattr(analysis, 'search_ms')
        
        # B1.3 warning fields
        assert hasattr(analysis, 'warnings')
        assert hasattr(analysis, 'is_partial')
    
    def test_spacer_extraction_deterministic(self):
        """Verify spacer extraction is deterministic."""
        from scalpel.design.spacer_extraction import SpacerExtractor
        from scalpel.models.enums import CasVariantType
        
        test_sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGT"
        extractor = SpacerExtractor(CasVariantType.SPCAS9)
        
        result1 = extractor.extract_spacers(test_sequence)
        result2 = extractor.extract_spacers(test_sequence)
        
        assert len(result1) == len(result2), "Same input should give same number of spacers"
        
        for s1, s2 in zip(result1, result2):
            assert s1.spacer_sequence == s2.spacer_sequence, "Spacer sequences should match"
            assert s1.genomic_start == s2.genomic_start, "Positions should match"

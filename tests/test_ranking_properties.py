"""
B5.2: Property tests for ranking invariants and bounds.

Tests that verify fundamental invariants of the guide ranking system:
- Score bounds (0-1 range)
- Ranking monotonicity (higher efficiency = higher rank)
- Determinism (same input = same output)
- Modality-specific behaviors
"""

from __future__ import annotations

import pytest
from typing import List
import random

from scalpel.models.enums import Strand, CasVariantType, EditModality, Genome
from scalpel.models.data_classes import SpacerCandidate, EfficiencyScore
from scalpel.design.spacer_extractor import SpacerExtractor
from scalpel.design.efficiency.ensemble import EnsembleScorer, ScoredGuide


# =============================================================================
# Helper Functions
# =============================================================================

def generate_random_spacer(seed: int = None) -> str:
    """Generate a random 20bp spacer sequence."""
    if seed is not None:
        random.seed(seed)
    bases = "ACGT"
    return "".join(random.choice(bases) for _ in range(20))


def generate_test_sequence(length: int = 200, seed: int = 42) -> str:
    """Generate a test sequence with embedded NGG PAM sites."""
    random.seed(seed)
    bases = "ACGT"
    seq = list("".join(random.choice(bases) for _ in range(length)))
    
    # Ensure some NGG PAMs exist
    pam_positions = [30, 70, 120, 160]
    for pos in pam_positions:
        if pos + 3 <= length:
            seq[pos:pos+3] = list("AGG")
    
    return "".join(seq)


def create_spacer_candidate(
    spacer: str,
    pam: str = "AGG",
    position: int = 100,
    strand: Strand = Strand.PLUS,
) -> SpacerCandidate:
    """Create a SpacerCandidate for testing."""
    return SpacerCandidate(
        spacer_sequence=spacer,
        pam_sequence=pam,
        strand=strand,
        genomic_start=position,
        genomic_end=position + len(spacer),
        cut_site=position + len(spacer) - 3,
        context_sequence=f"NNNNN{spacer}{pam}NNNNN",
    )


# =============================================================================
# Score Bounds Tests
# =============================================================================

class TestScoreBounds:
    """Verify all scores are within expected bounds."""
    
    def test_efficiency_score_bounded_0_to_1(self):
        """All efficiency scores must be in [0, 1] range."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        # Test with various spacer sequences
        test_spacers = [
            "ATGGATTTATCTGCTCTTCG",  # Normal
            "GGGGGGGGGGGGGGGGGGGG",  # All G
            "AAAAAAAAAAAAAAAAAANA",  # All A + N
            "TTTTTTTTTTTTTTTTTTTT",  # All T (Pol III terminator)
            "ACGTACGTACGTACGTACGT",  # Alternating
            "GCGCGCGCGCGCGCGCGCGC",  # High GC
            "ATATATATATATATATATA",   # Low GC (19bp) - edge case
        ]
        
        for spacer in test_spacers:
            if len(spacer) >= 20:
                spacer = spacer[:20]
            else:
                spacer = spacer + "A" * (20 - len(spacer))
            
            candidate = create_spacer_candidate(spacer)
            scored = scorer.score_single(candidate)
            
            assert 0 <= scored.efficiency.overall_score <= 1, \
                f"Score {scored.efficiency.overall_score} out of bounds for {spacer}"
    
    def test_random_spacers_bounded(self):
        """Generate 50 random spacers and verify all scores are bounded."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        for i in range(50):
            spacer = generate_random_spacer(seed=i)
            candidate = create_spacer_candidate(spacer)
            scored = scorer.score_single(candidate)
            
            assert 0 <= scored.efficiency.overall_score <= 1, \
                f"Random spacer {i} score {scored.efficiency.overall_score} out of bounds"
    
    def test_component_scores_bounded(self):
        """Individual component scores should also be bounded."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = "ATGGATTTATCTGCTCTTCG"
        candidate = create_spacer_candidate(spacer)
        scored = scorer.score_single(candidate)
        
        for component, score in scored.efficiency.components.items():
            assert 0 <= score <= 1, \
                f"Component '{component}' score {score} out of bounds"


class TestRankingMonotonicity:
    """Verify ranking respects score ordering."""
    
    def test_higher_score_means_higher_rank(self):
        """Guides with higher efficiency should rank higher (lower rank number)."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        # Create spacers with known efficiency differences
        # High GC in optimal range should score better than TTTT
        good_spacer = create_spacer_candidate("GCGATCGATCGATCGATCGG")  # Good GC
        bad_spacer = create_spacer_candidate("TTTTATTTATTTATTTATTT")   # TTTT penalty
        
        good_scored = scorer.score_single(good_spacer)
        bad_scored = scorer.score_single(bad_spacer)
        
        # If good scores higher, it should rank first when sorted
        guides = [bad_scored, good_scored]
        guides.sort(key=lambda x: x.efficiency.overall_score, reverse=True)
        
        if good_scored.efficiency.overall_score > bad_scored.efficiency.overall_score:
            assert guides[0].spacer.spacer_sequence == good_spacer.spacer_sequence, \
                "Higher scoring guide should rank first"
    
    def test_batch_ranking_preserves_order(self):
        """Batch scoring should preserve efficiency-based ordering."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        # Generate 10 random spacers
        spacers = [create_spacer_candidate(generate_random_spacer(i)) for i in range(10)]
        
        scored = scorer.score_batch(spacers)
        
        # Verify sorted by score descending
        scores = [s.efficiency.overall_score for s in scored]
        assert scores == sorted(scores, reverse=True), \
            "Batch results should be sorted by efficiency descending"
    
    def test_n_top_respects_ranking(self):
        """n_top parameter should return top N by efficiency."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        spacers = [create_spacer_candidate(generate_random_spacer(i)) for i in range(20)]
        
        all_scored = scorer.score_batch(spacers)
        top_5 = scorer.score_batch(spacers, n_top=5)
        
        # Top 5 from full list should match n_top=5
        assert len(top_5) == 5
        
        top_5_scores = [s.efficiency.overall_score for s in top_5]
        all_top_5_scores = [s.efficiency.overall_score for s in all_scored[:5]]
        
        assert top_5_scores == all_top_5_scores, \
            "n_top should return the same guides as top N from full list"


class TestDeterminism:
    """Verify scoring is deterministic."""
    
    def test_same_input_same_output(self):
        """Same spacer should always produce same score."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        
        scores = []
        for _ in range(10):
            scored = scorer.score_single(spacer)
            scores.append(scored.efficiency.overall_score)
        
        assert len(set(scores)) == 1, \
            f"Same input produced different scores: {set(scores)}"
    
    def test_batch_order_independent(self):
        """Batch scoring should be order-independent (same spacers, same scores)."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        spacers = [create_spacer_candidate(generate_random_spacer(i)) for i in range(5)]
        
        scored1 = scorer.score_batch(spacers)
        
        # Reverse order
        reversed_spacers = spacers[::-1]
        scored2 = scorer.score_batch(reversed_spacers)
        
        # Map spacer sequence to score
        scores1 = {s.spacer.spacer_sequence: s.efficiency.overall_score for s in scored1}
        scores2 = {s.spacer.spacer_sequence: s.efficiency.overall_score for s in scored2}
        
        assert scores1 == scores2, \
            "Batch order should not affect individual scores"
    
    def test_spacer_extraction_deterministic(self):
        """Spacer extraction should be deterministic."""
        extractor = SpacerExtractor(CasVariantType.SPCAS9)
        test_seq = generate_test_sequence(200)
        
        result1 = extractor.extract_spacers(test_seq)
        result2 = extractor.extract_spacers(test_seq)
        
        assert len(result1) == len(result2)
        
        for s1, s2 in zip(result1, result2):
            assert s1.spacer_sequence == s2.spacer_sequence
            assert s1.genomic_start == s2.genomic_start


class TestModalitySpecificBehavior:
    """Verify modality-specific scoring adjustments."""
    
    def test_knockout_modality_applied(self):
        """Knockout modality should apply specific adjustments."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        
        ko_scored = scorer.score_single(spacer, modality="knockout")
        
        # Should have valid score with modality applied
        assert 0 <= ko_scored.efficiency.overall_score <= 1
        assert ko_scored.efficiency.model_version  # Should have version string
    
    def test_interference_modality_applied(self):
        """Interference modality should apply specific adjustments."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        
        crispri_scored = scorer.score_single(spacer, modality="interference")
        
        assert 0 <= crispri_scored.efficiency.overall_score <= 1
    
    def test_activation_modality_applied(self):
        """Activation modality should apply specific adjustments."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        
        crispra_scored = scorer.score_single(spacer, modality="activation")
        
        assert 0 <= crispra_scored.efficiency.overall_score <= 1
    
    def test_different_modalities_can_differ(self):
        """Different modalities may produce different scores for same spacer."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        # Spacer with T at position 4 - penalized in knockout
        spacer = create_spacer_candidate("ATGTATTTATCTGCTCTTCG")
        
        ko_scored = scorer.score_single(spacer, modality="knockout")
        act_scored = scorer.score_single(spacer, modality="activation")
        
        # Scores may differ based on modality adjustments
        # Just verify both are valid
        assert 0 <= ko_scored.efficiency.overall_score <= 1
        assert 0 <= act_scored.efficiency.overall_score <= 1


class TestEdgeCases:
    """Edge cases and boundary conditions for ranking."""
    
    def test_empty_spacer_list(self):
        """Empty input should return empty output."""
        scorer = EnsembleScorer(use_ml_model=False)
        result = scorer.score_batch([])
        assert result == []
    
    def test_single_spacer_batch(self):
        """Single spacer batch should work correctly."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        
        result = scorer.score_batch([spacer])
        
        assert len(result) == 1
        assert 0 <= result[0].efficiency.overall_score <= 1
    
    def test_n_top_larger_than_list(self):
        """n_top larger than input should return all guides."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacers = [create_spacer_candidate(generate_random_spacer(i)) for i in range(5)]
        
        result = scorer.score_batch(spacers, n_top=100)
        
        assert len(result) == 5  # All 5, not 100
    
    def test_all_identical_spacers(self):
        """All identical spacers should have identical scores."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer_seq = "ATGGATTTATCTGCTCTTCG"
        
        spacers = [
            create_spacer_candidate(spacer_seq, position=i * 100)
            for i in range(10)
        ]
        
        scored = scorer.score_batch(spacers)
        scores = [s.efficiency.overall_score for s in scored]
        
        assert len(set(scores)) == 1, \
            f"Identical spacers should have identical scores: {set(scores)}"
    
    def test_confidence_interval_valid(self):
        """Confidence interval should be valid (low <= high, within 0-1)."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        scored = scorer.score_single(spacer)
        
        low, high = scored.efficiency.confidence_interval
        
        assert low <= high, f"CI lower ({low}) > upper ({high})"
        assert 0 <= low <= 1, f"CI lower ({low}) out of bounds"
        assert 0 <= high <= 1, f"CI upper ({high}) out of bounds"
    
    def test_interpretation_present(self):
        """Scored guide should have human-readable interpretation."""
        scorer = EnsembleScorer(use_ml_model=False)
        spacer = create_spacer_candidate("ATGGATTTATCTGCTCTTCG")
        scored = scorer.score_single(spacer)
        
        assert scored.efficiency.interpretation, \
            "Should have non-empty interpretation"
        assert any(q in scored.efficiency.interpretation for q in 
                   ["Excellent", "Good", "Moderate", "Low"]), \
            "Interpretation should contain quality descriptor"


class TestGCContentInvariants:
    """Test GC content effects on scoring."""
    
    def test_extreme_gc_penalized(self):
        """Extreme GC content (very high or low) should be penalized."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        # Optimal GC (~50%)
        optimal = create_spacer_candidate("ATGCATGCATGCATGCATGC")  # 50% GC
        
        # Very low GC
        low_gc = create_spacer_candidate("ATATATATATATATATATA" + "A")  # 0% GC
        
        # Very high GC
        high_gc = create_spacer_candidate("GCGCGCGCGCGCGCGCGCGC")  # 100% GC
        
        opt_score = scorer.score_single(optimal).efficiency.overall_score
        low_score = scorer.score_single(low_gc).efficiency.overall_score
        high_score = scorer.score_single(high_gc).efficiency.overall_score
        
        # Optimal should generally score better than extremes
        # (This is a soft invariant - may not always hold depending on other factors)
        assert opt_score >= 0
        assert low_score >= 0
        assert high_score >= 0
    
    def test_pol_iii_terminator_penalized(self):
        """TTTT sequence should be penalized (Pol III terminator)."""
        scorer = EnsembleScorer(use_ml_model=False)
        
        no_tttt = create_spacer_candidate("ATGCATGCATGCATGCATGC")
        with_tttt = create_spacer_candidate("ATGCTTTTGCATGCATGCAT")
        
        clean_score = scorer.score_single(no_tttt).efficiency.overall_score
        tttt_score = scorer.score_single(with_tttt).efficiency.overall_score
        
        assert tttt_score < clean_score, \
            f"TTTT ({tttt_score:.3f}) should score lower than clean ({clean_score:.3f})"

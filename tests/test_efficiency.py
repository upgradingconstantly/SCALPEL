"""Tests for efficiency scoring."""

from __future__ import annotations

import pytest

from scalpel.design.efficiency.ensemble import EnsembleScorer
from scalpel.design.efficiency.crispron_model import CRISPRonPredictor
from scalpel.design.efficiency.rule_based import RuleBasedScorer
from scalpel.models.data_classes import SpacerCandidate


class TestRuleBasedScoring:
    """Tests for rule-based efficiency scoring."""

    def test_gc_content_optimal_range(self) -> None:
        scorer = RuleBasedScorer()
        _, details = scorer.score("GCGCATATATGCGCATATAT", "AGG")
        assert details["gc_content"].score == 1.0

    def test_gc_content_low_penalty(self) -> None:
        scorer = RuleBasedScorer()
        _, details = scorer.score("ATATATATATATATATATAT", "AGG")
        assert details["gc_content"].score < 0.5

    def test_gc_content_high_penalty(self) -> None:
        scorer = RuleBasedScorer()
        _, details = scorer.score("GCGCGCGCGCGCGCGCGCGC", "AGG")
        assert details["gc_content"].score < 0.5

    def test_seed_gc_optimal(self) -> None:
        scorer = RuleBasedScorer()
        _, details = scorer.score("GCGCATATATATATATATAT", "AGG")
        assert details["gc_in_seed"].score == 1.0

    def test_homopolymer_tttt_penalty(self) -> None:
        scorer = RuleBasedScorer()
        _, details = scorer.score("AAAATTTTAAAATTTTAAAA", "AGG")
        assert details["tt_terminator"].score == 0.2

    def test_homopolymer_gggg_penalty(self) -> None:
        scorer = RuleBasedScorer()
        _, details = scorer.score("AAAAGGGGAAAAGGGGAAAA", "AGG")
        assert details["homopolymer"].score <= 0.5

    def test_position_20_preference(self) -> None:
        scorer = RuleBasedScorer()
        _, favored = scorer.score("ATCGATCGATCGATCGATCG", "AGG")
        _, unfavored = scorer.score("ATCGATCGATCGATCGATCT", "AGG")
        assert favored["position_20"].score > unfavored["position_20"].score

    def test_position_1_preference(self) -> None:
        scorer = RuleBasedScorer()
        _, favored = scorer.score("GTCGATCGATCGATCGATCG", "AGG")
        _, unfavored = scorer.score("ATCGATCGATCGATCGATCG", "AGG")
        assert favored["position_1"].score > unfavored["position_1"].score


class TestMLScoring:
    """Tests for CRISPRon predictor behavior."""

    def test_model_output_range(self) -> None:
        predictor = CRISPRonPredictor(use_rule_fallback=True)
        score = predictor.predict("ATCGATCGATCGATCGATCGAGG")
        assert 0.0 <= score <= 1.0

    def test_batch_prediction(self) -> None:
        predictor = CRISPRonPredictor(use_rule_fallback=True)
        scores = predictor.predict_batch(
            ["ATCGATCGATCGATCGATCGAGG", "GCGCGCGCGCGCGCGCGCGCAGG"]
        )
        assert len(scores) == 2
        assert all(0.0 <= s <= 1.0 for s in scores)

    def test_confidence_interval(self, sample_spacer_candidate: SpacerCandidate) -> None:
        scored = EnsembleScorer(use_ml_model=False).score_single(sample_spacer_candidate)
        lo, hi = scored.efficiency.confidence_interval
        assert 0.0 <= lo <= hi <= 1.0


class TestEnsembleScoring:
    """Tests for ensemble efficiency scorer."""

    def test_ensemble_combines_scores(self, sample_spacer_candidate: SpacerCandidate) -> None:
        class _FakePredictor:
            def predict(self, sequence: str) -> float:
                assert sequence
                return 0.9

        scorer = EnsembleScorer(use_ml_model=True, rule_weight=0.4, ml_weight=0.6)
        scorer._crispron = _FakePredictor()

        rule_only, _ = scorer.rule_scorer.score(
            sample_spacer_candidate.spacer_sequence,
            sample_spacer_candidate.pam_sequence,
        )
        scored = scorer.score_single(sample_spacer_candidate)
        expected = 0.4 * rule_only + 0.6 * 0.9
        assert scored.efficiency.overall_score == pytest.approx(expected, abs=1e-6)

    def test_score_bounds(self, sample_spacer_candidate: SpacerCandidate) -> None:
        scored = EnsembleScorer(use_ml_model=False).score_single(sample_spacer_candidate)
        assert 0.0 <= scored.efficiency.overall_score <= 1.0

    def test_interpretation_generated(self, sample_spacer_candidate: SpacerCandidate) -> None:
        scored = EnsembleScorer(use_ml_model=False).score_single(sample_spacer_candidate)
        assert scored.efficiency.interpretation

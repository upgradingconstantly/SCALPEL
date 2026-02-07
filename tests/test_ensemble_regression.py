"""Regression tests for ensemble efficiency scoring."""

from __future__ import annotations

import pytest

from scalpel.design.efficiency.ensemble import EnsembleScorer
from scalpel.design.efficiency.features import FeatureExtractor, SpacerFeatures
from scalpel.models.data_classes import SpacerCandidate


def test_score_single_no_longer_depends_on_missing_ml_model_attr(
    sample_spacer_candidate: SpacerCandidate,
) -> None:
    scorer = EnsembleScorer(use_ml_model=True)
    scorer._crispron = None
    assert not hasattr(scorer, "_ml_model")

    scored = scorer.score_single(sample_spacer_candidate)
    assert 0.0 <= scored.efficiency.overall_score <= 1.0
    assert scored.efficiency.model_version == "SCALPEL Ensemble v1 (rule)"


def test_predict_ml_fallback_uses_deterministic_position_bases() -> None:
    spacer = "ATCGATCGATCGATCGATCG"
    pam = "AGG"
    features = FeatureExtractor().extract(spacer, pam)

    scorer = EnsembleScorer(use_ml_model=True)
    scorer._crispron = None
    ml_fallback = scorer._predict_ml(features)
    expected_rule_score = scorer.rule_scorer.score(spacer, pam)[0]

    assert ml_fallback == pytest.approx(expected_rule_score, abs=1e-9)
    assert ml_fallback != 0.5


def test_predict_ml_empty_features_returns_neutral() -> None:
    features = SpacerFeatures(spacer="", pam="")
    scorer = EnsembleScorer(use_ml_model=False)
    assert scorer._predict_ml(features) == 0.5

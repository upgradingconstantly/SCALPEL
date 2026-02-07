"""
Ensemble efficiency scorer.

Combines rule-based and ML-based scoring for robust efficiency prediction.
"""

from __future__ import annotations

from typing import List, Optional, Dict
from dataclasses import dataclass

from scalpel.models.data_classes import SpacerCandidate, EfficiencyScore
from scalpel.design.efficiency.features import FeatureExtractor, SpacerFeatures
from scalpel.design.efficiency.rule_based import RuleBasedScorer, RuleScore
from scalpel.design.efficiency.crispron_model import CRISPRonPredictor, get_predictor


@dataclass
class ScoredGuide:
    """A spacer with its efficiency scores."""
    spacer: SpacerCandidate
    efficiency: EfficiencyScore
    features: SpacerFeatures
    rule_scores: Dict[str, RuleScore]


class EnsembleScorer:
    """
    Ensemble efficiency scorer combining multiple scoring methods.
    
    Methods:
    1. Rule-based scoring (Doench rules)
    2. Position-specific scoring
    3. Modality-specific adjustments
    4. CRISPRon ML model (CNN-based deep learning)
    """
    
    def __init__(
        self,
        use_ml_model: bool = True,
        rule_weight: float = 0.5,
        ml_weight: float = 0.5,
    ):
        self.use_ml_model = use_ml_model
        self.rule_weight = rule_weight
        self.ml_weight = ml_weight
        
        self.feature_extractor = FeatureExtractor()
        self.rule_scorer = RuleBasedScorer()
        
        # CRISPRon ML model (lazy loaded global instance)
        self._crispron: CRISPRonPredictor | None = None
        if self.use_ml_model:
            try:
                self._crispron = get_predictor()
            except Exception:
                self._crispron = None
    
    def score_single(
        self,
        spacer: SpacerCandidate,
        modality: str = "knockout",
    ) -> ScoredGuide:
        """
        Score a single spacer candidate.
        
        Args:
            spacer: SpacerCandidate to score
            modality: Editing modality for specific adjustments
        
        Returns:
            ScoredGuide with all scoring details
        """
        # Extract features
        features = self.feature_extractor.extract(
            spacer.spacer_sequence,
            spacer.pam_sequence,
            spacer.context_sequence,
        )
        
        # Get rule-based score
        rule_score, rule_details = self.rule_scorer.score(
            spacer.spacer_sequence,
            spacer.pam_sequence,
        )
        
        # Blend in ML score when a predictor is available.
        if self.use_ml_model and self._crispron is not None:
            ml_score = self._predict_ml(features)
            combined_score = (
                self.rule_weight * rule_score + 
                self.ml_weight * ml_score
            )
        else:
            combined_score = rule_score
        
        # Apply modality-specific adjustments
        adjusted_score = self._apply_modality_adjustment(
            combined_score, features, modality
        )
        
        # Generate interpretation
        interpretation = self._generate_interpretation(
            adjusted_score, features, rule_details
        )
        
        # Create efficiency score (matching data_classes.py schema)
        components = {
            name: rs.score for name, rs in rule_details.items()
        }
        
        model_version = "SCALPEL Ensemble v1 (rule)"
        if self.use_ml_model and self._crispron is not None:
            model_version = "SCALPEL Ensemble v1 (rule+ml)"

        efficiency = EfficiencyScore(
            overall_score=adjusted_score,
            components=components,
            confidence_interval=(max(0, adjusted_score - 0.1), min(1, adjusted_score + 0.1)),
            interpretation=interpretation,
            model_version=model_version,
        )
        
        return ScoredGuide(
            spacer=spacer,
            efficiency=efficiency,
            features=features,
            rule_scores=rule_details,
        )
    
    def score_batch(
        self,
        spacers: List[SpacerCandidate],
        modality: str = "knockout",
        n_top: Optional[int] = None,
    ) -> List[ScoredGuide]:
        """
        Score and rank multiple spacers.
        
        Args:
            spacers: List of SpacerCandidate objects
            modality: Editing modality
            n_top: If set, return only top N guides
        
        Returns:
            List of ScoredGuide objects, sorted by efficiency
        """
        scored = []
        for spacer in spacers:
            scored_guide = self.score_single(spacer, modality)
            scored.append(scored_guide)
        
        # Sort by efficiency score (descending)
        scored.sort(key=lambda x: x.efficiency.overall_score, reverse=True)
        
        if n_top is not None:
            scored = scored[:n_top]
        
        return scored
    
    def _predict_ml(self, features: SpacerFeatures) -> float:
        """ML model prediction using CRISPRon-style deep learning."""
        spacer = "".join(
            features.position_bases.get(i, "N")
            for i in range(1, len(features.position_bases) + 1)
        )
        if not spacer:
            return 0.5

        pam = features.pam if features.pam else "NGG"

        if self._crispron is None:
            # No model is available; return rule-based fallback.
            return self.rule_scorer.score(spacer, pam)[0]

        # CRISPRonPredictor accepts 20-30bp and pads/truncates internally.
        sequence_for_model = f"{spacer}{features.pam}" if features.pam else spacer
        return self._crispron.predict(sequence_for_model)
    
    def _apply_modality_adjustment(
        self,
        score: float,
        features: SpacerFeatures,
        modality: str,
    ) -> float:
        """
        Apply modality-specific score adjustments.
        
        Different modalities have different optimal guide properties:
        - Knockout: Prefer high cutting efficiency
        - CRISPRi: Position relative to TSS matters most
        - CRISPRa: Distance from TSS matters
        - Base editing: Editing window position matters
        """
        adjustment = 0.0
        
        if modality == "knockout":
            # For knockout, penalize T at position 4 more
            if features.position_bases.get(4) == "T":
                adjustment -= 0.05
        
        elif modality in ("interference", "activation"):
            # For CRISPRi/a, self-complementarity matters less
            if features.self_complementarity_score > 0.5:
                adjustment += 0.05  # Reduce penalty
        
        return max(0, min(1, score + adjustment))
    
    def _generate_interpretation(
        self,
        score: float,
        features: SpacerFeatures,
        rule_scores: Dict[str, RuleScore],
    ) -> str:
        """Generate human-readable interpretation of the score."""
        if score >= 0.8:
            quality = "Excellent"
        elif score >= 0.6:
            quality = "Good"
        elif score >= 0.4:
            quality = "Moderate"
        else:
            quality = "Low"
        
        # Find key factors
        factors = []
        
        # Check for notable issues
        if rule_scores.get("tt_terminator") and rule_scores["tt_terminator"].score < 0.5:
            factors.append("contains Pol III terminator")
        if rule_scores.get("homopolymer") and rule_scores["homopolymer"].score < 0.6:
            factors.append("homopolymer run")
        if features.gc_content < 0.35:
            factors.append("low GC")
        elif features.gc_content > 0.75:
            factors.append("high GC")
        
        # Build interpretation
        interpretation = f"{quality} efficiency (score: {score:.2f})"
        if factors:
            interpretation += f". Notable: {', '.join(factors)}"
        
        return interpretation


def score_guides(
    spacers: List[SpacerCandidate],
    modality: str = "knockout",
    n_top: Optional[int] = None,
) -> List[ScoredGuide]:
    """Convenience function to score and rank guides."""
    scorer = EnsembleScorer()
    return scorer.score_batch(spacers, modality, n_top)

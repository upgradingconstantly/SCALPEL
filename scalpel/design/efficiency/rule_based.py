"""
Rule-based efficiency scoring.

Implements scoring rules based on Doench 2014/2016 studies and other 
published guidelines for gRNA efficiency prediction.
"""

from __future__ import annotations

from typing import Dict, Tuple
from dataclasses import dataclass
import math

from scalpel.design.efficiency.features import SpacerFeatures, FeatureExtractor


@dataclass
class RuleScore:
    """Score from a single rule with explanation."""
    rule_name: str
    score: float  # 0-1
    weight: float
    reason: str


class RuleBasedScorer:
    """
    Rule-based gRNA efficiency scorer.
    
    Implements rules from:
    - Doench et al. 2014 (Rule Set 1)
    - Doench et al. 2016 (Rule Set 2)
    - Moreno-Mateos et al. 2015
    - Other published guidelines
    """
    
    # Default weights for each rule
    DEFAULT_WEIGHTS = {
        "gc_content": 0.15,
        "gc_in_seed": 0.10,
        "homopolymer": 0.15,
        "position_20": 0.10,
        "position_1": 0.05,
        "position_bias": 0.10,
        "self_comp": 0.10,
        "pam_strength": 0.10,
        "dinucleotide": 0.05,
        "tt_terminator": 0.10,
    }
    
    def __init__(self, weights: Dict[str, float] = None):
        self.weights = weights or self.DEFAULT_WEIGHTS
        self.feature_extractor = FeatureExtractor()
    
    def score(self, spacer: str, pam: str = "NGG") -> Tuple[float, Dict[str, RuleScore]]:
        """
        Score a spacer sequence.
        
        Args:
            spacer: 20bp spacer sequence
            pam: PAM sequence
        
        Returns:
            Tuple of (overall_score, detailed_scores)
        """
        features = self.feature_extractor.extract(spacer, pam)
        
        scores = {}
        
        # 1. GC content (optimal: 40-70%)
        scores["gc_content"] = self._score_gc_content(features)
        
        # 2. GC in seed region (optimal: 4-8 G/C in positions 1-12)
        scores["gc_in_seed"] = self._score_gc_seed(features)
        
        # 3. Homopolymer penalty
        scores["homopolymer"] = self._score_homopolymers(features)
        
        # 4. Position 20 (G or C preferred)
        scores["position_20"] = self._score_position_20(features)
        
        # 5. Position 1 (G preferred)
        scores["position_1"] = self._score_position_1(features)
        
        # 6. General position bias
        scores["position_bias"] = self._score_position_bias(features)
        
        # 7. Self-complementarity (secondary structure)
        scores["self_comp"] = self._score_self_complementarity(features)
        
        # 8. PAM strength
        scores["pam_strength"] = self._score_pam_strength(features)
        
        # 9. Dinucleotide bias
        scores["dinucleotide"] = self._score_dinucleotides(features)
        
        # 10. Pol III terminator (TTTT)
        scores["tt_terminator"] = self._score_tt_terminator(features)
        
        # Calculate weighted average
        total_weight = sum(self.weights.get(name, 0) for name in scores)
        if total_weight == 0:
            overall = 0.5
        else:
            overall = sum(
                scores[name].score * self.weights.get(name, 0)
                for name in scores
            ) / total_weight
        
        return overall, scores
    
    def _score_gc_content(self, features: SpacerFeatures) -> RuleScore:
        """GC content scoring - optimal is 40-70%."""
        gc = features.gc_content
        
        if 0.4 <= gc <= 0.7:
            score = 1.0
            reason = f"GC content {gc:.0%} is optimal (40-70%)"
        elif 0.3 <= gc < 0.4 or 0.7 < gc <= 0.8:
            # Gaussian falloff
            if gc < 0.4:
                score = 0.8 - (0.4 - gc) * 2
            else:
                score = 0.8 - (gc - 0.7) * 2
            reason = f"GC content {gc:.0%} slightly outside optimal range"
        else:
            # Very low or high GC
            score = max(0.2, 1 - abs(gc - 0.55) * 2)
            reason = f"GC content {gc:.0%} suboptimal"
        
        return RuleScore("gc_content", score, self.weights["gc_content"], reason)
    
    def _score_gc_seed(self, features: SpacerFeatures) -> RuleScore:
        """GC in seed region (positions 1-12)."""
        gc_count = int(features.gc_in_seed * 12)
        
        if 4 <= gc_count <= 8:
            score = 1.0
            reason = f"{gc_count} G/C in seed region (optimal: 4-8)"
        elif 3 <= gc_count <= 9:
            score = 0.8
            reason = f"{gc_count} G/C in seed region (acceptable)"
        else:
            score = 0.5
            reason = f"{gc_count} G/C in seed region (suboptimal)"
        
        return RuleScore("gc_in_seed", score, self.weights["gc_in_seed"], reason)
    
    def _score_homopolymers(self, features: SpacerFeatures) -> RuleScore:
        """Penalize homopolymer runs."""
        max_len = features.max_homopolymer_length
        
        if max_len <= 2:
            score = 1.0
            reason = "No significant homopolymers"
        elif max_len == 3:
            score = 0.8
            reason = "Contains 3bp homopolymer"
        elif max_len == 4:
            score = 0.5
            reason = "Contains 4bp homopolymer (may affect expression)"
        else:
            score = 0.2
            reason = f"Contains {max_len}bp homopolymer (likely problematic)"
        
        return RuleScore("homopolymer", score, self.weights["homopolymer"], reason)
    
    def _score_position_20(self, features: SpacerFeatures) -> RuleScore:
        """Position 20 preference (G or C)."""
        base = features.position_bases.get(20, "N")
        
        if base in ("G", "C"):
            score = 1.0
            reason = f"{base} at position 20 (favorable)"
        elif base == "A":
            score = 0.7
            reason = "A at position 20 (neutral)"
        else:  # T
            score = 0.5
            reason = "T at position 20 (less favorable)"
        
        return RuleScore("position_20", score, self.weights["position_20"], reason)
    
    def _score_position_1(self, features: SpacerFeatures) -> RuleScore:
        """Position 1 preference (G)."""
        base = features.position_bases.get(1, "N")
        
        if base == "G":
            score = 1.0
            reason = "G at position 1 (favorable for U6 expression)"
        elif base == "A":
            score = 0.8
            reason = "A at position 1 (acceptable)"
        else:
            score = 0.6
            reason = f"{base} at position 1 (may reduce expression)"
        
        return RuleScore("position_1", score, self.weights["position_1"], reason)
    
    def _score_position_bias(self, features: SpacerFeatures) -> RuleScore:
        """General position-specific base preferences."""
        spacer = features.spacer
        
        score = 1.0
        penalties = []
        
        # Unfavorable positions from Doench studies
        unfavorable = {
            4: "T",   # T at position 4
            14: "T",  # T at position 14
            17: "C",  # C at position 17
        }
        
        for pos, bad_base in unfavorable.items():
            if len(spacer) >= pos and spacer[pos - 1] == bad_base:
                score -= 0.1
                penalties.append(f"{bad_base}{pos}")
        
        score = max(0.5, score)
        
        if penalties:
            reason = f"Unfavorable bases at: {', '.join(penalties)}"
        else:
            reason = "No unfavorable position-specific bases"
        
        return RuleScore("position_bias", score, self.weights["position_bias"], reason)
    
    def _score_self_complementarity(self, features: SpacerFeatures) -> RuleScore:
        """Secondary structure penalty."""
        comp_score = features.self_complementarity_score
        
        if comp_score < 0.3:
            score = 1.0
            reason = "Low self-complementarity"
        elif comp_score < 0.5:
            score = 0.7
            reason = "Moderate self-complementarity"
        else:
            score = 0.4
            reason = "High self-complementarity (may form hairpin)"
        
        return RuleScore("self_comp", score, self.weights["self_comp"], reason)
    
    def _score_pam_strength(self, features: SpacerFeatures) -> RuleScore:
        """PAM sequence strength."""
        strength = features.pam_strength
        
        if strength >= 0.95:
            reason = f"Strong PAM ({features.pam})"
        elif strength >= 0.8:
            reason = f"Good PAM ({features.pam})"
        else:
            reason = f"Weak PAM ({features.pam})"
        
        return RuleScore("pam_strength", strength, self.weights["pam_strength"], reason)
    
    def _score_dinucleotides(self, features: SpacerFeatures) -> RuleScore:
        """Dinucleotide composition bias."""
        dinucs = features.dinucleotide_counts
        
        score = 1.0
        issues = []
        
        # TT dinucleotides are generally unfavorable
        tt_count = dinucs.get("TT", 0)
        if tt_count >= 3:
            score -= 0.2
            issues.append(f"{tt_count} TT")
        
        # GG can be favorable in some positions
        gg_count = dinucs.get("GG", 0)
        if gg_count >= 3:
            score -= 0.1
            issues.append(f"{gg_count} GG")
        
        score = max(0.5, score)
        
        if issues:
            reason = f"Dinucleotide issues: {', '.join(issues)}"
        else:
            reason = "Balanced dinucleotide composition"
        
        return RuleScore("dinucleotide", score, self.weights["dinucleotide"], reason)
    
    def _score_tt_terminator(self, features: SpacerFeatures) -> RuleScore:
        """Pol III terminator sequence (TTTT)."""
        if features.has_tttt:
            score = 0.2
            reason = "Contains TTTT (Pol III terminator - may truncate)"
        else:
            score = 1.0
            reason = "No Pol III terminator sequence"
        
        return RuleScore("tt_terminator", score, self.weights["tt_terminator"], reason)


def score_spacer(spacer: str, pam: str = "NGG") -> Tuple[float, Dict[str, RuleScore]]:
    """Convenience function to score a spacer."""
    return RuleBasedScorer().score(spacer, pam)

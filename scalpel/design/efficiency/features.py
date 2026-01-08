"""
Sequence feature extraction for efficiency scoring.

Extracts features from spacer sequences used by rule-based and ML models.
"""

from __future__ import annotations

from typing import Dict, List, Tuple
from dataclasses import dataclass, field
import math


@dataclass
class SpacerFeatures:
    """Extracted features from a spacer sequence."""
    
    # Sequence
    spacer: str
    pam: str
    context: str = ""
    
    # GC content features
    gc_content: float = 0.0
    gc_in_seed: float = 0.0  # Positions 1-12
    gc_in_nonseed: float = 0.0  # Positions 13-20
    
    # Position-specific bases
    position_bases: Dict[int, str] = field(default_factory=dict)
    
    # Homopolymer features  
    max_homopolymer_length: int = 0
    has_tttt: bool = False
    has_gggg: bool = False
    homopolymer_positions: List[Tuple[int, int, str]] = field(default_factory=list)
    
    # Dinucleotide features
    dinucleotide_counts: Dict[str, int] = field(default_factory=dict)
    
    # Thermodynamic features (simplified)
    predicted_tm: float = 0.0
    self_complementarity_score: float = 0.0
    
    # PAM features
    pam_strength: float = 1.0


class FeatureExtractor:
    """Extracts sequence features for efficiency prediction."""
    
    # Position weights from Doench 2014/2016 studies
    FAVORABLE_BASES = {
        1: "G",    # G at position 1 is favorable
        2: None,
        3: None,
        4: None,
        16: None,
        17: None,
        18: None,
        19: None,
        20: "G",   # G or C at position 20 is favorable
    }
    
    UNFAVORABLE_BASES = {
        4: "T",    # T at position 4 is unfavorable
        14: "T",   # T at position 14
        17: "C",   # C at position 17
    }
    
    def extract(self, spacer: str, pam: str = "", context: str = "") -> SpacerFeatures:
        """
        Extract all features from a spacer sequence.
        
        Args:
            spacer: 20bp spacer sequence
            pam: PAM sequence (e.g., "NGG")
            context: Extended sequence context
        
        Returns:
            SpacerFeatures dataclass with all computed features
        """
        spacer = spacer.upper()
        pam = pam.upper()
        
        features = SpacerFeatures(
            spacer=spacer,
            pam=pam,
            context=context,
        )
        
        # GC content
        features.gc_content = self._calc_gc_content(spacer)
        features.gc_in_seed = self._calc_gc_content(spacer[:12])
        features.gc_in_nonseed = self._calc_gc_content(spacer[12:])
        
        # Position-specific bases
        features.position_bases = {i + 1: base for i, base in enumerate(spacer)}
        
        # Homopolymers
        homopolymers = self._find_homopolymers(spacer)
        features.homopolymer_positions = homopolymers
        features.max_homopolymer_length = max((h[1] - h[0] for h in homopolymers), default=0)
        features.has_tttt = "TTTT" in spacer
        features.has_gggg = "GGGG" in spacer
        
        # Dinucleotides
        features.dinucleotide_counts = self._count_dinucleotides(spacer)
        
        # Thermodynamics
        features.predicted_tm = self._predict_tm(spacer)
        features.self_complementarity_score = self._calc_self_complementarity(spacer)
        
        # PAM strength
        features.pam_strength = self._score_pam(pam)
        
        return features
    
    def _calc_gc_content(self, seq: str) -> float:
        """Calculate GC content as fraction."""
        if len(seq) == 0:
            return 0.0
        gc_count = seq.count("G") + seq.count("C")
        return gc_count / len(seq)
    
    def _find_homopolymers(self, seq: str, min_length: int = 3) -> List[Tuple[int, int, str]]:
        """Find homopolymer runs (3+ identical bases)."""
        homopolymers = []
        i = 0
        while i < len(seq):
            base = seq[i]
            j = i + 1
            while j < len(seq) and seq[j] == base:
                j += 1
            length = j - i
            if length >= min_length:
                homopolymers.append((i, j, base))
            i = j
        return homopolymers
    
    def _count_dinucleotides(self, seq: str) -> Dict[str, int]:
        """Count all dinucleotide occurrences."""
        counts = {}
        for i in range(len(seq) - 1):
            dinuc = seq[i:i + 2]
            counts[dinuc] = counts.get(dinuc, 0) + 1
        return counts
    
    def _predict_tm(self, seq: str) -> float:
        """Simplified melting temperature prediction."""
        gc_count = seq.count("G") + seq.count("C")
        at_count = seq.count("A") + seq.count("T")
        
        if len(seq) < 14:
            # Wallace rule for short oligos
            return 2 * at_count + 4 * gc_count
        else:
            # Simplified nearest-neighbor approximation
            return 64.9 + 41 * (gc_count - 16.4) / len(seq)
    
    def _calc_self_complementarity(self, seq: str) -> float:
        """
        Calculate self-complementarity score.
        
        Higher score = more likely to form secondary structures.
        """
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        
        max_score = 0
        # Check for potential hairpin structures
        for i in range(len(seq) - 6):
            for j in range(i + 4, len(seq)):
                # Count complementary bases in potential stem
                matches = 0
                for k in range(min(i + 1, len(seq) - j)):
                    if i - k < 0 or j + k >= len(seq):
                        break
                    if seq[i - k] == complement.get(seq[j + k], ""):
                        matches += 1
                    else:
                        break
                max_score = max(max_score, matches)
        
        # Normalize to 0-1
        return min(1.0, max_score / 6)
    
    def _score_pam(self, pam: str) -> float:
        """Score PAM strength based on known preferences."""
        # For SpCas9, NGG is strongest
        pam_scores = {
            "AGG": 1.0,
            "TGG": 1.0,
            "CGG": 0.95,
            "GGG": 0.95,
            "NAG": 0.3,  # Weak PAM
        }
        
        for pattern, score in pam_scores.items():
            if self._pam_matches(pam, pattern):
                return score
        
        return 0.8  # Default for unknown PAMs
    
    def _pam_matches(self, pam: str, pattern: str) -> bool:
        """Check if PAM matches pattern (with N wildcards)."""
        if len(pam) != len(pattern):
            return False
        for p, pat in zip(pam, pattern):
            if pat != "N" and p != pat:
                return False
        return True


def extract_features(spacer: str, pam: str = "", context: str = "") -> SpacerFeatures:
    """Convenience function to extract features."""
    return FeatureExtractor().extract(spacer, pam, context)

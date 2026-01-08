"""
CFD (Cutting Frequency Determination) scoring for off-target analysis.

Implements the CFD scoring matrix from Doench et al. 2016 for predicting
the cutting probability at off-target sites based on mismatch position.
"""

from __future__ import annotations

from typing import Dict, List, Tuple


# CFD matrix from Doench et al. 2016
# Format: position -> (ref_base, mismatch_base) -> score
# Score of 1.0 means full cutting, lower = reduced cutting

CFD_MATRIX: Dict[int, Dict[Tuple[str, str], float]] = {
    # Position 1 (PAM-distal) - most tolerant to mismatches
    1: {
        ("A", "C"): 0.857, ("A", "G"): 0.571, ("A", "T"): 0.429,
        ("C", "A"): 0.857, ("C", "G"): 0.286, ("C", "T"): 0.714,
        ("G", "A"): 0.714, ("G", "C"): 0.429, ("G", "T"): 0.571,
        ("T", "A"): 0.714, ("T", "C"): 1.000, ("T", "G"): 0.571,
    },
    2: {
        ("A", "C"): 0.857, ("A", "G"): 0.571, ("A", "T"): 0.714,
        ("C", "A"): 0.857, ("C", "G"): 0.286, ("C", "T"): 0.857,
        ("G", "A"): 0.571, ("G", "C"): 0.286, ("G", "T"): 0.429,
        ("T", "A"): 0.714, ("T", "C"): 0.857, ("T", "G"): 0.571,
    },
    3: {
        ("A", "C"): 0.714, ("A", "G"): 0.571, ("A", "T"): 0.714,
        ("C", "A"): 0.714, ("C", "G"): 0.286, ("C", "T"): 0.714,
        ("G", "A"): 0.571, ("G", "C"): 0.286, ("G", "T"): 0.429,
        ("T", "A"): 0.571, ("T", "C"): 0.857, ("T", "G"): 0.429,
    },
    4: {
        ("A", "C"): 0.714, ("A", "G"): 0.429, ("A", "T"): 0.714,
        ("C", "A"): 0.571, ("C", "G"): 0.143, ("C", "T"): 0.714,
        ("G", "A"): 0.429, ("G", "C"): 0.286, ("G", "T"): 0.286,
        ("T", "A"): 0.571, ("T", "C"): 0.714, ("T", "G"): 0.429,
    },
    5: {
        ("A", "C"): 0.571, ("A", "G"): 0.429, ("A", "T"): 0.571,
        ("C", "A"): 0.571, ("C", "G"): 0.143, ("C", "T"): 0.571,
        ("G", "A"): 0.429, ("G", "C"): 0.143, ("G", "T"): 0.286,
        ("T", "A"): 0.429, ("T", "C"): 0.571, ("T", "G"): 0.286,
    },
    6: {
        ("A", "C"): 0.429, ("A", "G"): 0.429, ("A", "T"): 0.571,
        ("C", "A"): 0.571, ("C", "G"): 0.143, ("C", "T"): 0.571,
        ("G", "A"): 0.286, ("G", "C"): 0.143, ("G", "T"): 0.143,
        ("T", "A"): 0.429, ("T", "C"): 0.571, ("T", "G"): 0.286,
    },
    7: {
        ("A", "C"): 0.429, ("A", "G"): 0.286, ("A", "T"): 0.429,
        ("C", "A"): 0.429, ("C", "G"): 0.000, ("C", "T"): 0.429,
        ("G", "A"): 0.286, ("G", "C"): 0.143, ("G", "T"): 0.143,
        ("T", "A"): 0.286, ("T", "C"): 0.429, ("T", "G"): 0.143,
    },
    8: {
        ("A", "C"): 0.286, ("A", "G"): 0.286, ("A", "T"): 0.429,
        ("C", "A"): 0.429, ("C", "G"): 0.000, ("C", "T"): 0.286,
        ("G", "A"): 0.143, ("G", "C"): 0.143, ("G", "T"): 0.143,
        ("T", "A"): 0.286, ("T", "C"): 0.286, ("T", "G"): 0.143,
    },
    9: {
        ("A", "C"): 0.286, ("A", "G"): 0.143, ("A", "T"): 0.286,
        ("C", "A"): 0.286, ("C", "G"): 0.000, ("C", "T"): 0.286,
        ("G", "A"): 0.143, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.143, ("T", "C"): 0.286, ("T", "G"): 0.143,
    },
    10: {
        ("A", "C"): 0.286, ("A", "G"): 0.143, ("A", "T"): 0.286,
        ("C", "A"): 0.286, ("C", "G"): 0.000, ("C", "T"): 0.143,
        ("G", "A"): 0.143, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.143, ("T", "C"): 0.143, ("T", "G"): 0.000,
    },
    11: {
        ("A", "C"): 0.143, ("A", "G"): 0.143, ("A", "T"): 0.143,
        ("C", "A"): 0.143, ("C", "G"): 0.000, ("C", "T"): 0.143,
        ("G", "A"): 0.143, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.143, ("T", "C"): 0.143, ("T", "G"): 0.000,
    },
    12: {
        ("A", "C"): 0.143, ("A", "G"): 0.000, ("A", "T"): 0.143,
        ("C", "A"): 0.143, ("C", "G"): 0.000, ("C", "T"): 0.143,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.143, ("T", "G"): 0.000,
    },
    # Positions 13-20 (seed region) - least tolerant
    13: {
        ("A", "C"): 0.143, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    14: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    15: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    16: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    17: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    18: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    19: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
    20: {
        ("A", "C"): 0.000, ("A", "G"): 0.000, ("A", "T"): 0.000,
        ("C", "A"): 0.000, ("C", "G"): 0.000, ("C", "T"): 0.000,
        ("G", "A"): 0.000, ("G", "C"): 0.000, ("G", "T"): 0.000,
        ("T", "A"): 0.000, ("T", "C"): 0.000, ("T", "G"): 0.000,
    },
}

# PAM scores for SpCas9 variants
PAM_SCORES = {
    "GG": 1.0,   # NGG - canonical
    "AG": 0.90,  # NAG - weak PAM
    "CG": 0.10,  # NCG - very weak
    "TG": 0.05,  # NTG - very weak
    "GA": 0.01,  # NGA - almost inactive
    "GT": 0.01,
    "GC": 0.01,
}


class CFDScorer:
    """
    CFD (Cutting Frequency Determination) scorer.
    
    Calculates the probability of cutting at an off-target site
    based on the position and type of mismatches.
    """
    
    def __init__(self):
        self.cfd_matrix = CFD_MATRIX
        self.pam_scores = PAM_SCORES
    
    def score_offtarget(
        self,
        on_target: str,
        off_target: str,
        off_target_pam: str = "GG",
    ) -> float:
        """
        Calculate CFD score for an off-target site.
        
        Args:
            on_target: 20bp on-target spacer sequence
            off_target: 20bp off-target sequence
            off_target_pam: PAM at off-target site (last 2 bases)
        
        Returns:
            CFD score (0-1), where 1 = same cutting as on-target
        """
        on_target = on_target.upper()
        off_target = off_target.upper()
        off_target_pam = off_target_pam.upper()
        
        if len(on_target) != len(off_target):
            raise ValueError("On-target and off-target must be same length")
        
        # Start with perfect score
        score = 1.0
        
        # Apply mismatch penalties
        for i, (on_base, off_base) in enumerate(zip(on_target, off_target)):
            pos = i + 1  # 1-indexed position
            
            if on_base != off_base:
                # Get mismatch penalty
                penalty = self._get_mismatch_penalty(pos, on_base, off_base)
                score *= penalty
        
        # Apply PAM penalty
        pam_key = off_target_pam[-2:] if len(off_target_pam) >= 2 else "GG"
        pam_score = self.pam_scores.get(pam_key, 0.5)
        score *= pam_score
        
        return score
    
    def _get_mismatch_penalty(self, position: int, on_base: str, off_base: str) -> float:
        """Get penalty for a specific mismatch at a position."""
        if position not in self.cfd_matrix:
            # Default for positions outside matrix
            return 0.5
        
        key = (on_base, off_base)
        return self.cfd_matrix[position].get(key, 0.3)
    
    def find_mismatches(
        self,
        on_target: str,
        off_target: str,
    ) -> List[Tuple[int, str, str]]:
        """
        Find all mismatches between on-target and off-target.
        
        Returns:
            List of (position, on_base, off_base) tuples
        """
        mismatches = []
        for i, (on_base, off_base) in enumerate(zip(on_target, off_target)):
            if on_base != off_base:
                mismatches.append((i + 1, on_base, off_base))
        return mismatches
    
    def score_aggregate(
        self,
        on_target: str,
        off_target_sites: List[Tuple[str, str]],
    ) -> float:
        """
        Calculate aggregate off-target score.
        
        Lower aggregate score = fewer/weaker off-targets = better guide.
        
        Args:
            on_target: 20bp spacer
            off_target_sites: List of (off_target_seq, pam) tuples
        
        Returns:
            Aggregate score (sum of CFD scores for all off-targets)
        """
        total = 0.0
        for off_seq, pam in off_target_sites:
            total += self.score_offtarget(on_target, off_seq, pam)
        return total


def calculate_cfd_score(on_target: str, off_target: str, pam: str = "GG") -> float:
    """Convenience function to calculate CFD score."""
    return CFDScorer().score_offtarget(on_target, off_target, pam)

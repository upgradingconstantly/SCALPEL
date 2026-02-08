"""
Dual-Nickase Designer (Paired Cas9n).

Designs paired gRNAs for high-specificity editing using Cas9 nickase (D10A).
Two guides on opposite strands create offset nicks that form a DSB.
"""

from typing import Any, List, Tuple, Optional

from scalpel.core.plugin_registry import ModalityPlugin, register_modality
from scalpel.models.enums import EditModality, CasVariantType, Strand
from scalpel.models.data_classes import (
    ResolvedTarget,
    DesignedGuide,
    DesignResults,
    SpacerCandidate,
)


@register_modality(EditModality.DUAL_NICKASE)
class DualNickaseDesigner(ModalityPlugin):
    """
    Design paired gRNAs for dual-nickase editing.
    
    Strategy:
    1. Find pairs of guides on opposite strands
    2. Optimal spacing: 30-100 bp between nicks
    3. PAM-out orientation preferred (less off-target)
    4. D10A nickase creates 5' overhangs
    
    Benefits:
    - Dramatically reduced off-target effects
    - Only paired nicks create DSB
    - Single nicks are repaired by high-fidelity BER
    """
    
    display_name = "Dual-Nickase (Paired Cas9n)"
    description = "Design paired guides for high-specificity editing with Cas9 nickase"
    supported_cas_variants = [
        CasVariantType.SPCAS9,  # Used as D10A nickase
        CasVariantType.SPCAS9_NG,
        CasVariantType.SPCAS9_VQR,
    ]
    
    # Optimal spacing between nicks
    MIN_NICK_SPACING = 30
    MAX_NICK_SPACING = 100
    OPTIMAL_SPACING = (40, 70)
    
    scoring_weights = {
        "efficiency_pair": 0.35,  # Combined efficiency of pair
        "spacing": 0.30,          # Distance between nicks
        "orientation": 0.20,      # PAM-out preferred
        "specificity": 0.15,
    }
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,  # Returns n_guides PAIRS
        **kwargs: Any,
    ) -> DesignResults:
        """
        Design paired nickase guides.
        
        Returns guide pairs optimized for dual-nickase editing.
        """
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        
        cas_variant = target.specification.cas_variant
        
        # Extract all spacers
        extractor = SpacerExtractor(cas_variant)
        all_spacers = extractor.extract_spacers(
            target.sequence,
            chromosome=target.chromosome,
            start_position=target.start,
        )
        
        # Separate by strand
        plus_strand = [s for s in all_spacers if s.strand == Strand.PLUS]
        minus_strand = [s for s in all_spacers if s.strand == Strand.MINUS]
        
        # Find valid pairs
        valid_pairs = self._find_valid_pairs(plus_strand, minus_strand)
        
        # Score pairs
        scorer = EnsembleScorer()
        scored_pairs: List[Tuple[DesignedGuide, DesignedGuide, float]] = []
        
        for guide1, guide2, spacing in valid_pairs:
            # Score each guide individually
            scored1 = scorer.score_single(guide1, modality="knockout")
            scored2 = scorer.score_single(guide2, modality="knockout")
            
            # Combined efficiency (geometric mean)
            combined_eff = (scored1.efficiency.overall_score * 
                          scored2.efficiency.overall_score) ** 0.5
            
            # Spacing score
            spacing_score = self._score_spacing(spacing)
            
            # Orientation score (PAM-out is preferred)
            orientation_score = self._score_orientation(guide1, guide2)
            
            # Combined pair score
            pair_score = (
                self.scoring_weights["efficiency_pair"] * combined_eff +
                self.scoring_weights["spacing"] * spacing_score +
                self.scoring_weights["orientation"] * orientation_score +
                self.scoring_weights["specificity"] * 0.7
            )
            
            dg1 = DesignedGuide(
                spacer=guide1,
                efficiency_score=scored1.efficiency,
                composite_score=pair_score,
                position_score=spacing_score,
                rank=0,
            )
            dg2 = DesignedGuide(
                spacer=guide2,
                efficiency_score=scored2.efficiency,
                composite_score=pair_score,
                position_score=spacing_score,
                rank=0,
            )
            
            scored_pairs.append((dg1, dg2, pair_score))
        
        # Sort by pair score
        scored_pairs.sort(key=lambda x: x[2], reverse=True)
        
        # Return top pairs as alternating guides (guide1, guide2, guide1, guide2...)
        result_guides: List[DesignedGuide] = []
        for i, (g1, g2, score) in enumerate(scored_pairs[:n_guides], 1):
            g1.rank = i * 2 - 1  # Odd ranks
            g2.rank = i * 2      # Even ranks
            result_guides.append(g1)
            result_guides.append(g2)
        
        return DesignResults(
            target=target,
            modality=EditModality.DUAL_NICKASE,
            cas_variant=cas_variant,
            guides=result_guides[:n_guides * 2],  # Return pairs
            parameters={
                "n_pairs": n_guides, 
                "min_spacing": self.MIN_NICK_SPACING,
                "max_spacing": self.MAX_NICK_SPACING,
                **kwargs
            },
        )
    
    def _find_valid_pairs(
        self,
        plus_strand: List[SpacerCandidate],
        minus_strand: List[SpacerCandidate],
    ) -> List[Tuple[SpacerCandidate, SpacerCandidate, int]]:
        """Find guide pairs with valid nick spacing."""
        pairs = []
        
        for guide_plus in plus_strand:
            for guide_minus in minus_strand:
                # Calculate spacing between cut sites
                spacing = abs(guide_plus.cut_site - guide_minus.cut_site)
                
                if self.MIN_NICK_SPACING <= spacing <= self.MAX_NICK_SPACING:
                    pairs.append((guide_plus, guide_minus, spacing))
        
        return pairs
    
    def _score_spacing(self, spacing: int) -> float:
        """Score based on nick spacing (optimal: 40-70 bp)."""
        opt_min, opt_max = self.OPTIMAL_SPACING
        
        if opt_min <= spacing <= opt_max:
            return 1.0
        elif spacing < opt_min:
            return 0.7 + 0.3 * (spacing - self.MIN_NICK_SPACING) / (opt_min - self.MIN_NICK_SPACING)
        else:  # spacing > opt_max
            return 0.7 + 0.3 * (self.MAX_NICK_SPACING - spacing) / (self.MAX_NICK_SPACING - opt_max)
    
    def _score_orientation(
        self,
        guide1: SpacerCandidate,
        guide2: SpacerCandidate,
    ) -> float:
        """
        Score PAM orientation.
        
        PAM-out (5'-N20-PAM...PAM-N20-5'): preferred, fewer off-targets
        PAM-in: acceptable but less specific
        """
        # For SpCas9, PAM is 3' of spacer
        # PAM-out: plus strand guide upstream, minus strand guide downstream
        # This creates 5' overhangs with D10A nickase
        
        if guide1.strand == Strand.PLUS:
            plus_guide, minus_guide = guide1, guide2
        else:
            plus_guide, minus_guide = guide2, guide1
        
        # PAM-out: plus strand cut site < minus strand cut site
        if plus_guide.cut_site < minus_guide.cut_site:
            return 1.0  # PAM-out orientation
        else:
            return 0.7  # PAM-in orientation (still valid)
    
    def get_default_parameters(self) -> dict:
        return {
            "min_spacing": self.MIN_NICK_SPACING,
            "max_spacing": self.MAX_NICK_SPACING,
            "optimal_spacing": self.OPTIMAL_SPACING,
        }

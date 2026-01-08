"""
CRISPRa (CRISPR Activation) Designer.

Designs gRNAs for transcriptional activation using dCas9-VP64/VPR/SAM.
"""

from typing import Any, List

from scalpel.core.plugin_registry import ModalityPlugin, register_modality
from scalpel.models.enums import EditModality, CasVariantType, Strand
from scalpel.models.data_classes import (
    ResolvedTarget,
    DesignedGuide,
    DesignResults,
    SpacerCandidate,
)


@register_modality(EditModality.ACTIVATION)
class ActivationDesigner(ModalityPlugin):
    """
    Design gRNAs for CRISPRa (transcriptional activation).
    
    Strategy:
    1. Target -400 to -50 upstream of TSS (promoter region)
    2. Multiple guides for synergistic activation (SAM/SunTag)
    3. Avoid core promoter elements
    """
    
    display_name = "CRISPRa (Activation)"
    description = "Design guides for transcriptional activation using dCas9-VP64/VPR"
    supported_cas_variants = [
        CasVariantType.SPCAS9,  # Used as dCas9
        CasVariantType.SPCAS9_NG,
    ]
    
    # Optimal window upstream of TSS
    PROMOTER_START = -400
    PROMOTER_END = -50
    OPTIMAL_START = -200
    OPTIMAL_END = -100
    
    scoring_weights = {
        "efficiency": 0.30,
        "promoter_position": 0.45,  # Position in promoter is critical
        "strand": 0.10,
        "specificity": 0.15,
    }
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        **kwargs: Any,
    ) -> DesignResults:
        """
        Design CRISPRa guides targeting the promoter region.
        """
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        
        cas_variant = target.specification.cas_variant
        
        # Get TSS position
        tss = self._get_tss(target)
        
        # Extract spacers
        extractor = SpacerExtractor(cas_variant)
        all_spacers = extractor.extract_spacers(
            target.sequence,
            chromosome=target.chromosome,
            start_position=target.start,
        )
        
        # Filter to promoter window
        valid_spacers = self._filter_by_promoter(all_spacers, tss, target)
        
        # Score
        scorer = EnsembleScorer()
        scored_guides: List[DesignedGuide] = []
        
        for spacer in valid_spacers:
            scored = scorer.score_single(spacer, modality="activation")
            
            # Calculate promoter position score
            position_score = self._score_promoter_position(spacer, tss, target)
            
            # Strand score
            strand_score = 0.8  # Less critical for CRISPRa
            
            # Combined score
            composite = (
                self.scoring_weights["efficiency"] * scored.efficiency.overall_score +
                self.scoring_weights["promoter_position"] * position_score +
                self.scoring_weights["strand"] * strand_score +
                self.scoring_weights["specificity"] * 0.7
            )
            
            guide = DesignedGuide(
                spacer=spacer,
                efficiency_score=scored.efficiency,
                composite_score=composite,
                position_score=position_score,
                domain_score=strand_score,
                rank=0,
            )
            scored_guides.append(guide)
        
        # Sort and rank
        scored_guides.sort(key=lambda g: g.composite_score, reverse=True)
        for i, guide in enumerate(scored_guides[:n_guides], 1):
            guide.rank = i
        
        return DesignResults(
            target=target,
            modality=EditModality.ACTIVATION,
            cas_variant=cas_variant,
            guides=scored_guides[:n_guides],
            parameters={"n_guides": n_guides, "tss": tss, **kwargs},
        )
    
    def _get_tss(self, target: ResolvedTarget) -> int:
        """Get transcription start site."""
        if target.transcript:
            return target.transcript.tss
        return target.start if target.strand == Strand.PLUS else target.end
    
    def _filter_by_promoter(
        self,
        spacers: List[SpacerCandidate],
        tss: int,
        target: ResolvedTarget,
    ) -> List[SpacerCandidate]:
        """Filter spacers to promoter region."""
        valid = []
        for spacer in spacers:
            distance = spacer.cut_site - tss
            if target.strand == Strand.MINUS:
                distance = -distance
            
            if self.PROMOTER_START <= distance <= self.PROMOTER_END:
                valid.append(spacer)
        
        # If no promoter guides, use all (may be targeting enhancer)
        return valid if valid else spacers[:50]
    
    def _score_promoter_position(
        self,
        spacer: SpacerCandidate,
        tss: int,
        target: ResolvedTarget,
    ) -> float:
        """Score based on position in promoter."""
        distance = spacer.cut_site - tss
        if target.strand == Strand.MINUS:
            distance = -distance
        
        # Optimal range (-200 to -100)
        if self.OPTIMAL_START <= distance <= self.OPTIMAL_END:
            return 1.0
        # Acceptable range
        elif self.PROMOTER_START <= distance <= self.PROMOTER_END:
            return 0.7
        else:
            return 0.3
    
    def get_default_parameters(self) -> dict:
        return {
            "tss_upstream_min": 50,
            "tss_upstream_max": 400,
        }

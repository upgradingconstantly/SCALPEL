"""
CRISPRi (CRISPR Interference) Designer.

Designs gRNAs for transcriptional repression using dCas9-KRAB.
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


@register_modality(EditModality.INTERFERENCE)
class InterferenceDesigner(ModalityPlugin):
    """
    Design gRNAs for CRISPRi (transcriptional repression).
    
    Strategy:
    1. Target -50 to +300 relative to TSS (optimal window)
    2. Prefer template (non-coding) strand
    3. dCas9-KRAB fusions block transcription elongation
    """
    
    display_name = "CRISPRi (Interference)"
    description = "Design guides for transcriptional repression using dCas9-KRAB"
    supported_cas_variants = [
        CasVariantType.SPCAS9,  # Used as dCas9
        CasVariantType.SPCAS9_NG,
    ]
    
    # Optimal window relative to TSS
    TSS_WINDOW_START = -50
    TSS_WINDOW_END = 300
    OPTIMAL_START = 50
    OPTIMAL_END = 150
    
    scoring_weights = {
        "efficiency": 0.35,
        "tss_distance": 0.40,  # Distance from TSS is critical
        "strand": 0.15,       # Template strand preferred
        "specificity": 0.10,
    }
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        **kwargs: Any,
    ) -> DesignResults:
        """
        Design CRISPRi guides targeting near the TSS.
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
        
        # Filter to TSS window
        valid_spacers = self._filter_by_tss_window(all_spacers, tss, target)
        
        # Score efficiency
        scorer = EnsembleScorer()
        scored_guides: List[DesignedGuide] = []
        
        for spacer in valid_spacers:
            scored = scorer.score_single(spacer, modality="interference")
            
            # Calculate TSS distance score
            tss_score = self._score_tss_distance(spacer, tss)
            
            # Calculate strand score (template strand preferred)
            strand_score = self._score_strand(spacer, target)
            
            # Combined score
            composite = (
                self.scoring_weights["efficiency"] * scored.efficiency.overall_score +
                self.scoring_weights["tss_distance"] * tss_score +
                self.scoring_weights["strand"] * strand_score +
                self.scoring_weights["specificity"] * 0.7
            )
            
            guide = DesignedGuide(
                spacer=spacer,
                efficiency_score=scored.efficiency,
                composite_score=composite,
                position_score=tss_score,
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
            modality=EditModality.INTERFERENCE,
            cas_variant=cas_variant,
            guides=scored_guides[:n_guides],
            parameters={"n_guides": n_guides, "tss": tss, **kwargs},
        )
    
    def _get_tss(self, target: ResolvedTarget) -> int:
        """Get transcription start site."""
        if target.transcript:
            return target.transcript.tss
        return target.start if target.strand == Strand.PLUS else target.end
    
    def _filter_by_tss_window(
        self,
        spacers: List[SpacerCandidate],
        tss: int,
        target: ResolvedTarget,
    ) -> List[SpacerCandidate]:
        """Filter spacers to TSS window."""
        valid = []
        for spacer in spacers:
            distance = spacer.cut_site - tss
            if target.strand == Strand.MINUS:
                distance = -distance
            
            if self.TSS_WINDOW_START <= distance <= self.TSS_WINDOW_END:
                valid.append(spacer)
        
        return valid if valid else spacers[:50]  # Fallback
    
    def _score_tss_distance(self, spacer: SpacerCandidate, tss: int) -> float:
        """Score based on distance from TSS (optimal: +50 to +150)."""
        distance = abs(spacer.cut_site - tss)
        
        # Optimal range
        if self.OPTIMAL_START <= distance <= self.OPTIMAL_END:
            return 1.0
        # Acceptable range
        elif distance <= self.TSS_WINDOW_END:
            # Linear falloff
            if distance < self.OPTIMAL_START:
                return 0.7 + 0.3 * (distance / self.OPTIMAL_START)
            else:
                return 0.7 * (1 - (distance - self.OPTIMAL_END) / 
                             (self.TSS_WINDOW_END - self.OPTIMAL_END))
        else:
            return 0.3
    
    def _score_strand(
        self,
        spacer: SpacerCandidate,
        target: ResolvedTarget,
    ) -> float:
        """Template strand (non-coding) is preferred for CRISPRi."""
        if target.strand == Strand.PLUS:
            return 1.0 if spacer.strand == Strand.MINUS else 0.7
        else:
            return 1.0 if spacer.strand == Strand.PLUS else 0.7
    
    def get_default_parameters(self) -> dict:
        return {
            "tss_window_start": self.TSS_WINDOW_START,
            "tss_window_end": self.TSS_WINDOW_END,
            "prefer_template_strand": True,
        }

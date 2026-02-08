"""
CRISPRoff Designer (Epigenetic Silencing).

Designs gRNAs for heritable gene silencing via DNA methylation.
Uses dCas9-DNMT3A/DNMT3L fusions to deposit CpG methylation near TSS.
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


@register_modality(EditModality.CRISPROFF)
class CRISPRoffDesigner(ModalityPlugin):
    """
    Design gRNAs for CRISPRoff (epigenetic silencing).
    
    Strategy:
    1. Target CpG islands near TSS (-1000 to +500 bp)
    2. Prioritize CpG-rich regions for stable methylation
    3. dCas9-DNMT3A/DNMT3L deposits DNA methylation
    4. Heritable silencing without DNA sequence change
    """
    
    display_name = "CRISPRoff (Epigenetic Silencing)"
    description = "Design guides for heritable gene silencing via DNA methylation"
    supported_cas_variants = [
        CasVariantType.SPCAS9,  # Used as dCas9
        CasVariantType.SPCAS9_NG,
    ]
    
    # Target window relative to TSS
    TSS_WINDOW_START = -1000
    TSS_WINDOW_END = 500
    OPTIMAL_START = -500
    OPTIMAL_END = 200
    
    # CpG scoring
    MIN_CPG_DENSITY = 0.6  # Minimum CpG observed/expected ratio
    
    scoring_weights = {
        "efficiency": 0.25,
        "tss_distance": 0.30,
        "cpg_density": 0.35,  # CpG density is critical for stable silencing
        "specificity": 0.10,
    }
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        **kwargs: Any,
    ) -> DesignResults:
        """
        Design CRISPRoff guides targeting CpG islands near TSS.
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
        
        # Filter to TSS/CpG island window
        valid_spacers = self._filter_by_tss_window(all_spacers, tss, target)
        
        # Score efficiency
        scorer = EnsembleScorer()
        scored_guides: List[DesignedGuide] = []
        
        for spacer in valid_spacers:
            scored = scorer.score_single(spacer, modality="interference")  # Similar to CRISPRi
            
            # Calculate TSS distance score
            tss_score = self._score_tss_distance(spacer, tss, target)
            
            # Calculate CpG density score
            cpg_score = self._score_cpg_density(spacer)
            
            # Combined score
            composite = (
                self.scoring_weights["efficiency"] * scored.efficiency.overall_score +
                self.scoring_weights["tss_distance"] * tss_score +
                self.scoring_weights["cpg_density"] * cpg_score +
                self.scoring_weights["specificity"] * 0.7
            )
            
            guide = DesignedGuide(
                spacer=spacer,
                efficiency_score=scored.efficiency,
                composite_score=composite,
                position_score=tss_score,
                domain_score=cpg_score,
                rank=0,
            )
            scored_guides.append(guide)
        
        # Sort and rank
        scored_guides.sort(key=lambda g: g.composite_score, reverse=True)
        for i, guide in enumerate(scored_guides[:n_guides], 1):
            guide.rank = i
        
        return DesignResults(
            target=target,
            modality=EditModality.CRISPROFF,
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
        """Filter spacers to CpG island window near TSS."""
        valid = []
        for spacer in spacers:
            distance = spacer.cut_site - tss
            if target.strand == Strand.MINUS:
                distance = -distance
            
            if self.TSS_WINDOW_START <= distance <= self.TSS_WINDOW_END:
                valid.append(spacer)
        
        return valid if valid else spacers[:50]  # Fallback
    
    def _score_tss_distance(
        self,
        spacer: SpacerCandidate,
        tss: int,
        target: ResolvedTarget,
    ) -> float:
        """Score based on distance from TSS."""
        distance = spacer.cut_site - tss
        if target.strand == Strand.MINUS:
            distance = -distance
        
        # Optimal range (-500 to +200)
        if self.OPTIMAL_START <= distance <= self.OPTIMAL_END:
            return 1.0
        # Acceptable range
        elif self.TSS_WINDOW_START <= distance <= self.TSS_WINDOW_END:
            if distance < self.OPTIMAL_START:
                return 0.6 + 0.4 * (1 - abs(distance - self.OPTIMAL_START) / 
                                    abs(self.TSS_WINDOW_START - self.OPTIMAL_START))
            else:
                return 0.6 + 0.4 * (1 - abs(distance - self.OPTIMAL_END) / 
                                    abs(self.TSS_WINDOW_END - self.OPTIMAL_END))
        else:
            return 0.3
    
    def _score_cpg_density(self, spacer: SpacerCandidate) -> float:
        """
        Score based on CpG density around the spacer.
        
        CpG density = (CpG count / expected CpG) 
        where expected = (C count * G count) / sequence length
        """
        # Use context sequence if available, otherwise spacer + PAM
        seq = spacer.context_sequence if spacer.context_sequence else spacer.spacer_sequence
        seq = seq.upper()
        
        if len(seq) < 10:
            return 0.5
        
        # Count CpG dinucleotides
        cpg_count = seq.count("CG")
        c_count = seq.count("C")
        g_count = seq.count("G")
        
        # Calculate observed/expected ratio
        expected = (c_count * g_count) / len(seq) if len(seq) > 0 else 0
        
        if expected > 0:
            observed_expected_ratio = cpg_count / expected
        else:
            observed_expected_ratio = 0
        
        # Score: higher CpG density = better for stable methylation
        if observed_expected_ratio >= 0.8:  # CpG island criterion
            return 1.0
        elif observed_expected_ratio >= 0.6:
            return 0.8
        elif observed_expected_ratio >= 0.4:
            return 0.6
        else:
            return 0.4
    
    def get_default_parameters(self) -> dict:
        return {
            "tss_window_start": self.TSS_WINDOW_START,
            "tss_window_end": self.TSS_WINDOW_END,
            "min_cpg_density": self.MIN_CPG_DENSITY,
        }

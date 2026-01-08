"""
CRISPR Knockout Designer.

Designs gRNAs for gene knockout via NHEJ-mediated frameshift.
"""

from typing import Any, List, Optional

from scalpel.core.plugin_registry import ModalityPlugin, register_modality
from scalpel.models.enums import EditModality, CasVariantType, Strand
from scalpel.models.data_classes import (
    ResolvedTarget,
    DesignedGuide,
    DesignResults,
    SpacerCandidate,
    EfficiencyScore,
)


@register_modality(EditModality.KNOCKOUT)
class KnockoutDesigner(ModalityPlugin):
    """
    Design gRNAs for gene knockout.
    
    Strategy:
    1. Target early, constitutive exons
    2. Prefer functional domains
    3. Avoid NMD-escape regions (last exon, last 50bp of penultimate exon)
    """
    
    display_name = "CRISPR Knockout"
    description = "Design guides for gene knockout via NHEJ-mediated frameshift"
    supported_cas_variants = [
        CasVariantType.SPCAS9,
        CasVariantType.SPCAS9_NG,
        CasVariantType.SACAS9,
        CasVariantType.CAS12A,
    ]
    
    # Scoring weights optimized for knockout
    scoring_weights = {
        "efficiency": 0.40,
        "position": 0.30,  # Early in CDS preferred
        "domain": 0.15,    # Targeting functional domains
        "specificity": 0.15,
    }
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        **kwargs: Any,
    ) -> DesignResults:
        """
        Design knockout guides for the given target.
        
        Full pipeline:
        1. Extract spacers from target region
        2. Filter by valid exons (skip first/last)
        3. Score efficiency
        4. Score position (early in CDS = better)
        5. Combine scores and rank
        """
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        
        cas_variant = target.specification.cas_variant
        params = {**self.get_default_parameters(), **kwargs}
        
        # Extract spacers
        extractor = SpacerExtractor(cas_variant)
        all_spacers = extractor.extract_spacers(
            target.sequence,
            chromosome=target.chromosome,
            start_position=target.start,
        )
        
        # Filter spacers by position (skip first/last exons)
        valid_spacers = self._filter_by_exon(all_spacers, target, params)
        
        # Score efficiency
        scorer = EnsembleScorer()
        scored_guides: List[DesignedGuide] = []
        
        for spacer in valid_spacers:
            # Get efficiency score
            scored = scorer.score_single(spacer, modality="knockout")
            
            # Calculate position score
            position_score = self._score_position(spacer, target)
            
            # Calculate domain score (if targeting functional domain)
            domain_score = self._score_domain(spacer, target)
            
            # Combined score with weights
            composite = (
                self.scoring_weights["efficiency"] * scored.efficiency.overall_score +
                self.scoring_weights["position"] * position_score +
                self.scoring_weights["domain"] * domain_score +
                self.scoring_weights["specificity"] * 0.7  # Placeholder until off-target run
            )
            
            guide = DesignedGuide(
                spacer=spacer,
                efficiency_score=scored.efficiency,
                composite_score=composite,
                position_score=position_score,
                domain_score=domain_score,
                rank=0,  # Will be set after sorting
            )
            scored_guides.append(guide)
        
        # Sort by composite score and assign ranks
        scored_guides.sort(key=lambda g: g.composite_score, reverse=True)
        for i, guide in enumerate(scored_guides[:n_guides], 1):
            guide.rank = i
        
        return DesignResults(
            target=target,
            modality=EditModality.KNOCKOUT,
            cas_variant=cas_variant,
            guides=scored_guides[:n_guides],
            parameters={"n_guides": n_guides, **params},
        )
    
    def _filter_by_exon(
        self,
        spacers: List[SpacerCandidate],
        target: ResolvedTarget,
        params: dict,
    ) -> List[SpacerCandidate]:
        """Filter spacers to valid exons only."""
        if not target.transcript or not target.transcript.exons:
            return spacers
        
        exons = target.transcript.exons
        valid_spacers = []
        
        for spacer in spacers:
            cut_site = spacer.cut_site
            
            # Find which exon contains the cut site
            for i, exon in enumerate(exons):
                if exon.start <= cut_site <= exon.end:
                    # Skip first exon if configured
                    if params.get("skip_first_exon") and i == 0:
                        break
                    # Skip last exon if configured
                    if params.get("skip_last_exon") and i == len(exons) - 1:
                        break
                    # Check minimum exon size
                    if exon.end - exon.start < params.get("min_exon_size", 50):
                        break
                    
                    spacer.in_exon = True
                    spacer.exon_number = exon.exon_number
                    valid_spacers.append(spacer)
                    break
        
        # If no valid exons found, return all spacers
        return valid_spacers if valid_spacers else spacers
    
    def _score_position(
        self,
        spacer: SpacerCandidate,
        target: ResolvedTarget,
    ) -> float:
        """
        Score guide position for knockout effectiveness.
        
        Earlier in CDS = better (more protein truncation).
        """
        if not target.transcript or not target.transcript.cds_start:
            return 0.5
        
        cds_start = target.transcript.cds_start
        cds_end = target.transcript.cds_end or target.transcript.end
        cut_site = spacer.cut_site
        
        # Position in CDS (0 = start, 1 = end)
        cds_length = cds_end - cds_start
        if cds_length <= 0:
            return 0.5
        
        relative_position = abs(cut_site - cds_start) / cds_length
        relative_position = max(0, min(1, relative_position))
        
        # Prefer first 30% of CDS
        if relative_position < 0.3:
            return 1.0
        elif relative_position < 0.5:
            return 0.8
        elif relative_position < 0.7:
            return 0.6
        else:
            return 0.4
    
    def _score_domain(
        self,
        spacer: SpacerCandidate,
        target: ResolvedTarget,
    ) -> float:
        """Score based on whether cut site is in a functional domain."""
        if not target.gene_info:
            return 0.5
        
        # Get domains from annotation database
        from scalpel.genome.annotations import AnnotationDatabase
        
        try:
            db = AnnotationDatabase(target.specification.genome)
            domains = db.get_protein_domains(target.gene_info.gene_id)
            
            for domain in domains:
                if domain.genomic_start <= spacer.cut_site <= domain.genomic_end:
                    if domain.is_catalytic:
                        return 1.0  # Catalytic domains are best
                    elif domain.is_essential:
                        return 0.9
                    else:
                        return 0.7
        except Exception:
            pass
        
        return 0.5  # No domain targeting
    
    def score_position(
        self,
        guide: DesignedGuide,
        target: ResolvedTarget,
    ) -> float:
        """Public interface for position scoring."""
        return self._score_position(guide.spacer, target)
    
    def get_default_parameters(self) -> dict:
        return {
            "skip_first_exon": True,
            "skip_last_exon": True,
            "prefer_constitutive": True,
            "min_exon_size": 50,
        }


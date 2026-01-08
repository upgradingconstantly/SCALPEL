"""
Prime Editor Designer.

Designs pegRNAs for prime editing with PBS and RTT components.
"""

from typing import Any, Optional, List

from scalpel.core.plugin_registry import ModalityPlugin, register_modality
from scalpel.models.enums import EditModality, CasVariantType
from scalpel.models.data_classes import (
    ResolvedTarget,
    DesignedGuide,
    DesignResults,
    SpacerCandidate,
)


@register_modality(EditModality.PRIME_EDIT)
class PrimeEditorDesigner(ModalityPlugin):
    """
    Design pegRNAs for prime editing.
    
    Prime editing requires:
    1. Spacer sequence (targets nick site)
    2. PBS (primer binding site) - 8-15 nt
    3. RTT (reverse transcription template) - encodes edit
    4. Optional: nicking gRNA for PE3/PE3b
    """
    
    display_name = "Prime Editor"
    description = "Design pegRNAs for precise editing without DSBs"
    supported_cas_variants = [
        CasVariantType.SPCAS9,  # Uses Cas9 nickase (H840A)
    ]
    
    # Prime editing scoring weights
    scoring_weights = {
        "spacer_efficiency": 0.25,
        "pbs_quality": 0.25,
        "rtt_quality": 0.25,
        "structure_score": 0.15,
        "specificity": 0.10,
    }
    
    # PBS parameters
    pbs_length_range = (8, 15)
    optimal_pbs_tm = (50, 60)  # Celsius
    
    # RTT parameters
    rtt_length_range = (10, 30)
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        desired_edit: Optional[str] = None,
        **kwargs: Any,
    ) -> DesignResults:
        """Design pegRNAs for the given target and desired edit."""
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        
        cas_variant = target.specification.cas_variant
        params = {**self.get_default_parameters(), **kwargs}
        
        # Extract all spacers from target region
        extractor = SpacerExtractor(cas_variant)
        all_spacers = extractor.extract_spacers(
            target.sequence,
            chromosome=target.chromosome,
            start_position=target.start,
        )
        
        # Score efficiency and generate pegRNA components
        scorer = EnsembleScorer()
        scored_guides: List[DesignedGuide] = []
        
        for spacer in all_spacers:
            # Get base efficiency score
            scored = scorer.score_single(spacer, modality="knockout")  # Use base scoring
            
            # Design PBS variants
            nick_position = len(spacer.spacer_sequence)  # Simplified nick pos
            pbs_variants = self.design_pbs_variants(
                spacer.spacer_sequence,
                target.sequence,
                nick_position,
            )
            
            # Use best PBS score
            best_pbs = pbs_variants[0] if pbs_variants else {"score": 0.5, "sequence": "", "tm": 0}
            pbs_score = best_pbs.get("score", 0.5)
            
            # Score nick position relative to edit site
            position_score = self._score_nick_position(spacer.cut_site, target)
            
            # RTT quality (simplified - would depend on desired_edit in production)
            rtt_score = 0.7  # Placeholder
            
            # Combined score with weights
            composite = (
                self.scoring_weights["spacer_efficiency"] * scored.efficiency.overall_score +
                self.scoring_weights["pbs_quality"] * pbs_score +
                self.scoring_weights["rtt_quality"] * rtt_score +
                self.scoring_weights["structure_score"] * 0.7 +  # Placeholder
                self.scoring_weights["specificity"] * 0.7  # Placeholder
            )
            
            guide = DesignedGuide(
                spacer=spacer,
                efficiency_score=scored.efficiency,
                composite_score=composite,
                position_score=position_score,
                domain_score=pbs_score,  # Reuse for PBS quality
                rank=0,
            )
            scored_guides.append(guide)
        
        # Sort by composite score and assign ranks
        scored_guides.sort(key=lambda g: g.composite_score, reverse=True)
        for i, guide in enumerate(scored_guides[:n_guides], 1):
            guide.rank = i
        
        return DesignResults(
            target=target,
            modality=EditModality.PRIME_EDIT,
            cas_variant=cas_variant,
            guides=scored_guides[:n_guides],
            parameters={"n_guides": n_guides, "desired_edit": desired_edit, **params},
        )
    
    def _score_nick_position(self, cut_site: int, target: ResolvedTarget) -> float:
        """Score nick position relative to target center."""
        # Optimal: nick is central to target region
        target_center = (target.start + target.end) / 2
        distance = abs(cut_site - target_center)
        max_distance = (target.end - target.start) / 2
        
        if max_distance <= 0:
            return 0.5
        
        # Closer to center = better
        normalized = 1.0 - min(distance / max_distance, 1.0)
        return normalized
    
    def score_position(
        self,
        guide: DesignedGuide,
        target: ResolvedTarget,
    ) -> float:
        """Score nick position relative to edit site."""
        # Optimal: nick site 0-30bp from edit
        return 0.5  # Placeholder
    
    def design_pbs_variants(
        self,
        spacer: str,
        target_sequence: str,
        nick_position: int,
    ) -> List[dict]:
        """
        Design PBS (primer binding site) variants.
        
        The PBS binds to the target strand 3' of the nick.
        """
        pbs_variants = []
        
        for length in range(self.pbs_length_range[0], self.pbs_length_range[1] + 1):
            # PBS template is sequence downstream of nick
            pbs_template = target_sequence[nick_position:nick_position + length]
            pbs_seq = self._reverse_complement(pbs_template)
            
            # Calculate melting temperature
            tm = self._calculate_tm(pbs_seq)
            
            # Score
            if self.optimal_pbs_tm[0] <= tm <= self.optimal_pbs_tm[1]:
                tm_score = 1.0
            else:
                tm_score = max(0, 1 - abs(tm - 55) / 20)
            
            pbs_variants.append({
                "sequence": pbs_seq,
                "length": length,
                "tm": tm,
                "score": tm_score,
            })
        
        return sorted(pbs_variants, key=lambda x: x["score"], reverse=True)
    
    def _reverse_complement(self, seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(complement.get(base, "N") for base in reversed(seq.upper()))
    
    def _calculate_tm(self, seq: str) -> float:
        """Calculate melting temperature using nearest-neighbor method (simplified)."""
        # Simplified Tm calculation
        gc_count = seq.upper().count("G") + seq.upper().count("C")
        at_count = seq.upper().count("A") + seq.upper().count("T")
        
        if len(seq) < 14:
            # Wallace rule for short oligos
            return 2 * at_count + 4 * gc_count
        else:
            # Simplified nearest-neighbor approximation
            return 64.9 + 41 * (gc_count - 16.4) / len(seq)
    
    def get_default_parameters(self) -> dict:
        return {
            "pbs_length_min": 8,
            "pbs_length_max": 15,
            "rtt_length_min": 10,
            "rtt_length_max": 30,
            "design_nicking_guide": True,
            "pe3_distance_min": 40,
            "pe3_distance_max": 120,
        }

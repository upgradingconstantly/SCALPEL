"""
Base Editor Designer (CBE/ABE).

Designs gRNAs for base editing:
- CBE (Cytosine Base Editor): C→T conversion
- ABE (Adenine Base Editor): A→G conversion
"""

from typing import Any, List, Tuple, Optional

from scalpel.core.plugin_registry import ModalityPlugin, register_modality
from scalpel.models.enums import EditModality, CasVariantType
from scalpel.models.data_classes import (
    ResolvedTarget,
    DesignedGuide,
    DesignResults,
    SpacerCandidate,
)


class BaseEditorDesigner(ModalityPlugin):
    """
    Design gRNAs for base editing.
    
    Key constraints:
    - Edit must occur in editing window (CBE: pos 4-8, ABE: pos 4-7)
    - Minimize bystander edits (other C/A in window)
    - Consider codon context for amino acid change
    """
    
    supported_cas_variants = [
        CasVariantType.SPCAS9,
        CasVariantType.SPCAS9_NG,
    ]
    
    def __init__(self, editor_type: str = "CBE"):
        """
        Args:
            editor_type: "CBE" for cytosine, "ABE" for adenine
        """
        self.editor_type = editor_type
        
        if editor_type == "CBE":
            self.target_base = "C"
            self.product_base = "T"
            self.editing_window = (4, 8)
            self.window_weights = {4: 0.6, 5: 0.9, 6: 1.0, 7: 0.9, 8: 0.5}
        else:  # ABE
            self.target_base = "A"
            self.product_base = "G"
            self.editing_window = (4, 7)
            self.window_weights = {4: 0.5, 5: 0.9, 6: 1.0, 7: 0.7}
    
    @property
    def display_name(self) -> str:
        return f"Base Editor ({self.editor_type})"
    
    @property
    def description(self) -> str:
        return f"Design guides for {self.target_base}→{self.product_base} base editing"
    
    # Base editing scoring weights
    scoring_weights = {
        "efficiency": 0.30,
        "window_position": 0.30,
        "bystander_penalty": 0.25,
        "specificity": 0.15,
    }
    
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        desired_edit: Optional[str] = None,  # e.g., "c.123C>T"
        **kwargs: Any,
    ) -> DesignResults:
        """Design base editing guides for the given target."""
        from scalpel.design import SpacerExtractor
        from scalpel.design.efficiency import EnsembleScorer
        
        cas_variant = target.specification.cas_variant
        params = {**self.get_default_parameters(), **kwargs}
        max_bystanders = params.get("max_bystanders", 2)
        
        # Extract all spacers from target region
        extractor = SpacerExtractor(cas_variant)
        all_spacers = extractor.extract_spacers(
            target.sequence,
            chromosome=target.chromosome,
            start_position=target.start,
        )
        
        # Filter spacers where target base falls in editing window
        valid_spacers = self._filter_by_window(all_spacers, max_bystanders)
        
        # Score efficiency
        scorer = EnsembleScorer()
        scored_guides: List[DesignedGuide] = []
        
        modality_key = "base_cbe" if self.editor_type == "CBE" else "base_abe"
        
        for spacer in valid_spacers:
            # Get efficiency score
            scored = scorer.score_single(spacer, modality=modality_key)
            
            # Calculate window position score (how well target is positioned)
            window_score = self._score_window_position(spacer.spacer_sequence)
            
            # Count bystanders and calculate penalty
            bystander_count, bystander_positions = self.count_bystanders(spacer.spacer_sequence)
            bystander_penalty = 1.0 - (bystander_count * 0.15)  # -15% per bystander
            bystander_penalty = max(0.2, bystander_penalty)  # Floor at 0.2
            
            # Combined score with weights
            composite = (
                self.scoring_weights["efficiency"] * scored.efficiency.overall_score +
                self.scoring_weights["window_position"] * window_score +
                self.scoring_weights["bystander_penalty"] * bystander_penalty +
                self.scoring_weights["specificity"] * 0.7  # Placeholder until off-target
            )
            
            guide = DesignedGuide(
                spacer=spacer,
                efficiency_score=scored.efficiency,
                composite_score=composite,
                position_score=window_score,
                domain_score=bystander_penalty,  # Reuse domain_score for bystander
                rank=0,
            )
            scored_guides.append(guide)
        
        # Sort by composite score and assign ranks
        scored_guides.sort(key=lambda g: g.composite_score, reverse=True)
        for i, guide in enumerate(scored_guides[:n_guides], 1):
            guide.rank = i
        
        return DesignResults(
            target=target,
            modality=EditModality.BASE_EDIT_CBE if self.editor_type == "CBE" else EditModality.BASE_EDIT_ABE,
            cas_variant=cas_variant,
            guides=scored_guides[:n_guides],
            parameters={"n_guides": n_guides, "desired_edit": desired_edit, **params},
        )
    
    def _filter_by_window(
        self,
        spacers: List[SpacerCandidate],
        max_bystanders: int,
    ) -> List[SpacerCandidate]:
        """Filter spacers to only those with target base in editing window."""
        valid = []
        window_start, window_end = self.editing_window
        
        for spacer in spacers:
            seq = spacer.spacer_sequence
            # Check if target base exists in editing window
            window_seq = seq[window_start - 1:window_end]  # Convert 1-indexed to 0-indexed
            
            if self.target_base in window_seq:
                # Count how many target bases (potential bystanders)
                target_count = window_seq.count(self.target_base)
                if target_count <= max_bystanders + 1:  # +1 for intended edit
                    valid.append(spacer)
        
        return valid
    
    def _score_window_position(self, spacer: str) -> float:
        """Score based on where target base falls in editing window."""
        window_start, window_end = self.editing_window
        window_seq = spacer[window_start - 1:window_end]
        
        # Find best position for target base
        best_score = 0.0
        for i, base in enumerate(window_seq):
            if base == self.target_base:
                pos = window_start + i
                weight = self.window_weights.get(pos, 0.5)
                best_score = max(best_score, weight)
        
        return best_score
    
    def score_position(
        self,
        guide: DesignedGuide,
        target: ResolvedTarget,
    ) -> float:
        """Score based on edit position in window and bystander count."""
        # This would be implemented with specific edit position
        return 0.5
    
    def count_bystanders(self, spacer: str) -> Tuple[int, List[int]]:
        """
        Count bystander bases in editing window.
        
        Returns:
            Tuple of (count, positions)
        """
        window_start, window_end = self.editing_window
        window_seq = spacer[window_start - 1:window_end]  # 1-indexed to 0-indexed
        
        positions = [
            i + window_start
            for i, base in enumerate(window_seq)
            if base == self.target_base
        ]
        
        return len(positions), positions
    
    def get_default_parameters(self) -> dict:
        return {
            "editor_type": self.editor_type,
            "max_bystanders": 2,
        }


@register_modality(EditModality.BASE_EDIT_CBE)
class CBEDesigner(BaseEditorDesigner):
    """Cytosine Base Editor (C→T) designer."""
    
    def __init__(self):
        super().__init__(editor_type="CBE")


@register_modality(EditModality.BASE_EDIT_ABE)
class ABEDesigner(BaseEditorDesigner):
    """Adenine Base Editor (A→G) designer."""
    
    def __init__(self):
        super().__init__(editor_type="ABE")

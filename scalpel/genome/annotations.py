"""
Genomic annotation database interface.

Provides lookup for regulatory elements, protein domains, and other annotations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass
from enum import Enum

from scalpel.models.enums import Genome


class AnnotationType(str, Enum):
    """Types of genomic annotations."""
    PROMOTER = "promoter"
    ENHANCER = "enhancer"
    SILENCER = "silencer"
    INSULATOR = "insulator"
    OPEN_CHROMATIN = "open_chromatin"
    PROTEIN_DOMAIN = "protein_domain"
    REPEAT = "repeat"
    CONSERVATION = "conservation"


@dataclass
class GenomicAnnotation:
    """A genomic annotation."""
    annotation_type: AnnotationType
    chromosome: str
    start: int
    end: int
    name: Optional[str] = None
    score: Optional[float] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def overlaps(self, start: int, end: int) -> bool:
        """Check if this annotation overlaps with a region."""
        return self.start < end and start < self.end


@dataclass
class ProteinDomain:
    """Protein domain annotation."""
    domain_id: str
    name: str
    gene_id: str
    start_aa: int  # Amino acid position
    end_aa: int
    genomic_start: int
    genomic_end: int
    is_catalytic: bool = False
    is_essential: bool = False
    is_binding: bool = False
    source: str = "Pfam"


@dataclass 
class RepeatElement:
    """Repeat element annotation."""
    repeat_class: str  # SINE, LINE, LTR, DNA, Simple_repeat, etc.
    repeat_name: str   # Alu, L1, etc.
    chromosome: str
    start: int
    end: int
    strand: str
    overlap_percent: float = 0.0


class AnnotationDatabase:
    """
    Interface for genomic annotations.
    
    Integrates with:
    - ENCODE regulatory elements
    - Pfam protein domains
    - RepeatMasker annotations
    - PhyloP conservation scores
    """
    
    def __init__(self, genome: Genome, data_dir: Optional[Path] = None):
        self.genome = genome
        
        if data_dir is None:
            from scalpel.config import get_config
            config = get_config()
            data_dir = config.data_dir / "annotations" / genome.value
        
        self.data_dir = data_dir
        self._loaded = False
    
    def get_annotations_at_position(
        self,
        chromosome: str,
        start: int,
        end: int,
        annotation_types: Optional[List[AnnotationType]] = None,
    ) -> List[GenomicAnnotation]:
        """
        Get annotations overlapping a genomic region.
        
        Args:
            chromosome: Chromosome name
            start: Start position
            end: End position
            annotation_types: Filter by annotation type (None = all)
        
        Returns:
            List of overlapping annotations
        """
        # TODO: Implement with actual annotation database
        # For now, return empty list
        return []
    
    def get_regulatory_elements(
        self,
        chromosome: str,
        position: int,
        window: int = 1000,
    ) -> List[GenomicAnnotation]:
        """Get regulatory elements near a position."""
        return self.get_annotations_at_position(
            chromosome,
            position - window,
            position + window,
            annotation_types=[
                AnnotationType.PROMOTER,
                AnnotationType.ENHANCER,
                AnnotationType.SILENCER,
            ]
        )
    
    def get_protein_domains(self, gene_id: str) -> List[ProteinDomain]:
        """
        Get protein domains for a gene.
        
        Args:
            gene_id: Ensembl gene ID
        
        Returns:
            List of protein domains
        """
        # TODO: Implement with Pfam/InterPro data
        # Demo data for testing
        if gene_id == "ENSG00000141510":  # TP53
            return [
                ProteinDomain(
                    domain_id="PF00870",
                    name="P53 DNA-binding domain",
                    gene_id=gene_id,
                    start_aa=94,
                    end_aa=292,
                    genomic_start=7674180,
                    genomic_end=7676272,
                    is_essential=True,
                    is_binding=True,
                ),
                ProteinDomain(
                    domain_id="PF07710",
                    name="P53 tetramerization domain",
                    gene_id=gene_id,
                    start_aa=323,
                    end_aa=356,
                    genomic_start=7673534,
                    genomic_end=7673836,
                    is_essential=True,
                ),
            ]
        return []
    
    def get_repeat_elements(
        self,
        chromosome: str,
        start: int,
        end: int,
    ) -> List[RepeatElement]:
        """
        Get repeat elements overlapping a region.
        
        Args:
            chromosome: Chromosome name
            start: Start position
            end: End position
        
        Returns:
            List of repeat elements
        """
        # TODO: Implement with RepeatMasker data
        return []
    
    def get_conservation_score(
        self,
        chromosome: str,
        position: int,
    ) -> Optional[float]:
        """
        Get conservation score (PhyloP) at a position.
        
        Returns:
            Conservation score (-20 to +10) or None if not available
        """
        # TODO: Implement with PhyloP bigWig
        return None
    
    def is_essential_gene(self, gene_id: str) -> bool:
        """
        Check if a gene is essential.
        
        Uses DepMap essentiality data.
        """
        # TODO: Implement with DepMap data
        # For now, hardcode known essential genes
        essential_genes = {
            "ENSG00000141510",  # TP53
            "ENSG00000012048",  # BRCA1
            "ENSG00000183765",  # POLR2A
            "ENSG00000105968",  # RPS19
        }
        return gene_id in essential_genes
    
    def is_cancer_gene(self, gene_id: str) -> bool:
        """
        Check if a gene is a known cancer gene.
        
        Uses COSMIC Cancer Gene Census.
        """
        # Demo data
        cancer_genes = {
            "ENSG00000141510",  # TP53 - tumor suppressor
            "ENSG00000012048",  # BRCA1 - tumor suppressor
            "ENSG00000146648",  # EGFR - oncogene
            "ENSG00000133703",  # KRAS - oncogene
            "ENSG00000181555",  # SETD2 - tumor suppressor
        }
        return gene_id in cancer_genes

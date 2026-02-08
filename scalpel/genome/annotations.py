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
    - DepMap essentiality data
    """

    def __init__(self, genome: Genome, data_dir: Optional[Path] = None):
        self.genome = genome

        if data_dir is None:
            from scalpel.config import get_config
            config = get_config()
            data_dir = config.data_dir / "annotations" / genome.value

        self.data_dir = data_dir
        self._loaded = False
        self._repeat_cache: Dict[str, List[RepeatElement]] = {}
        self._conservation_cache: Dict[str, Dict[int, float]] = {}
        self._depmap_essentiality: Dict[str, float] = {}
        self._load_demo_data()

    def _load_demo_data(self) -> None:
        """Load demonstration data for common use cases."""
        # DepMap-style essentiality scores (lower = more essential, <-0.5 = essential)
        self._depmap_essentiality = {
            "ENSG00000141510": -0.89,  # TP53
            "ENSG00000012048": -0.76,  # BRCA1
            "ENSG00000139618": -0.72,  # BRCA2
            "ENSG00000183765": -1.45,  # POLR2A (very essential)
            "ENSG00000105968": -1.12,  # RPS19
            "ENSG00000146648": -0.23,  # EGFR (less essential)
            "ENSG00000133703": -0.31,  # KRAS
            "ENSG00000171862": -0.95,  # PTEN
            "ENSG00000134982": -0.88,  # APC
            "ENSG00000136997": -0.67,  # MYC
        }

        # Demo repeat elements for TP53 region (chr17)
        self._repeat_cache["chr17"] = [
            RepeatElement(
                repeat_class="SINE",
                repeat_name="AluSx",
                chromosome="chr17",
                start=7668000,
                end=7668300,
                strand="+",
            ),
            RepeatElement(
                repeat_class="LINE",
                repeat_name="L1MA4",
                chromosome="chr17",
                start=7690000,
                end=7691500,
                strand="-",
            ),
            RepeatElement(
                repeat_class="Simple_repeat",
                repeat_name="(CA)n",
                chromosome="chr17",
                start=7675500,
                end=7675600,
                strand="+",
            ),
        ]

        # Demo conservation scores for TP53 DNA-binding domain
        import random
        rng = random.Random(42)
        conservation = {}
        # Highly conserved DNA-binding domain region
        for pos in range(7674180, 7676272):
            conservation[pos] = rng.uniform(4.0, 9.5)  # High conservation
        # Less conserved intron regions
        for pos in range(7668421, 7674180):
            conservation[pos] = rng.uniform(-2.0, 3.0)
        self._conservation_cache["chr17"] = conservation

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
            List of repeat elements with overlap percentage calculated
        """
        repeats = self._repeat_cache.get(chromosome, [])
        overlapping = []

        for repeat in repeats:
            if repeat.start < end and start < repeat.end:
                # Calculate overlap percentage
                overlap_start = max(start, repeat.start)
                overlap_end = min(end, repeat.end)
                overlap_length = overlap_end - overlap_start
                query_length = end - start
                overlap_pct = (overlap_length / query_length) * 100 if query_length > 0 else 0

                overlapping.append(RepeatElement(
                    repeat_class=repeat.repeat_class,
                    repeat_name=repeat.repeat_name,
                    chromosome=repeat.chromosome,
                    start=repeat.start,
                    end=repeat.end,
                    strand=repeat.strand,
                    overlap_percent=round(overlap_pct, 1),
                ))

        return overlapping
    
    def get_conservation_score(
        self,
        chromosome: str,
        position: int,
    ) -> Optional[float]:
        """
        Get conservation score (PhyloP) at a position.

        Returns:
            Conservation score (-20 to +10) or None if not available.
            Higher positive values = more conserved.
        """
        chr_cache = self._conservation_cache.get(chromosome)
        if chr_cache:
            return chr_cache.get(position)
        return None

    def get_average_conservation(
        self,
        chromosome: str,
        start: int,
        end: int,
    ) -> Optional[float]:
        """
        Get average conservation score for a region.

        Returns:
            Average PhyloP score or None if no data available.
        """
        chr_cache = self._conservation_cache.get(chromosome)
        if not chr_cache:
            return None

        scores = [chr_cache[pos] for pos in range(start, end) if pos in chr_cache]
        if scores:
            return sum(scores) / len(scores)
        return None
    
    def is_essential_gene(self, gene_id: str) -> bool:
        """
        Check if a gene is essential.

        Uses DepMap essentiality data where score < -0.5 indicates essentiality.
        """
        score = self._depmap_essentiality.get(gene_id)
        if score is not None:
            return score < -0.5
        return False

    def get_essentiality_score(self, gene_id: str) -> Optional[float]:
        """
        Get the DepMap essentiality score for a gene.

        Returns:
            Essentiality score (lower = more essential, <-0.5 = essential)
            or None if not available.
        """
        return self._depmap_essentiality.get(gene_id)

    def get_essentiality_interpretation(self, gene_id: str) -> str:
        """
        Get human-readable essentiality interpretation.
        """
        score = self._depmap_essentiality.get(gene_id)
        if score is None:
            return "Unknown"
        if score < -1.0:
            return "Strongly Essential"
        if score < -0.5:
            return "Essential"
        if score < 0:
            return "Moderately Important"
        return "Non-Essential"
    
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

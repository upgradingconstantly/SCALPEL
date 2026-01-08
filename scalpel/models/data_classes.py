"""
Pydantic data classes for SCALPEL.

All core data structures used throughout the platform.
"""

from __future__ import annotations

from typing import Optional, Dict, List, Tuple, Any
from pydantic import BaseModel, Field, computed_field
from datetime import datetime

from scalpel.models.enums import (
    TargetType,
    Genome,
    EditModality,
    CasVariantType,
    PAMPosition,
    Strand,
    RedFlagSeverity,
)


# =============================================================================
# Genomic Data Structures
# =============================================================================

class GenomicRegion(BaseModel):
    """A region in the genome."""
    chromosome: str
    start: int
    end: int
    strand: Strand = Strand.PLUS
    
    @computed_field
    @property
    def length(self) -> int:
        return self.end - self.start
    
    def overlaps(self, other: "GenomicRegion") -> bool:
        """Check if this region overlaps with another."""
        if self.chromosome != other.chromosome:
            return False
        return self.start < other.end and other.start < self.end
    
    def contains(self, position: int) -> bool:
        """Check if position is within this region."""
        return self.start <= position < self.end


class Exon(BaseModel):
    """Exon annotation."""
    exon_number: int
    start: int
    end: int
    is_constitutive: bool = True
    
    @computed_field
    @property
    def length(self) -> int:
        return self.end - self.start


class Transcript(BaseModel):
    """Transcript annotation."""
    transcript_id: str
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: Strand
    exons: List[Exon] = Field(default_factory=list)
    cds_start: Optional[int] = None
    cds_end: Optional[int] = None
    biotype: str = "protein_coding"
    is_mane_select: bool = False
    
    @computed_field
    @property
    def tss(self) -> int:
        """Transcription start site."""
        return self.start if self.strand == Strand.PLUS else self.end
    
    @computed_field
    @property
    def cds_length(self) -> int:
        """Length of coding sequence."""
        if self.cds_start is None or self.cds_end is None:
            return 0
        return self.cds_end - self.cds_start


class GeneInfo(BaseModel):
    """Gene annotation."""
    gene_id: str
    symbol: str
    chromosome: str
    start: int
    end: int
    strand: Strand
    biotype: str = "protein_coding"
    description: Optional[str] = None
    aliases: List[str] = Field(default_factory=list)
    transcripts: List[Transcript] = Field(default_factory=list)


# =============================================================================
# Target Specification
# =============================================================================

class TargetSpecification(BaseModel):
    """Complete specification for a CRISPR target."""
    target_type: TargetType
    target_value: str
    genome: Genome
    modality: EditModality
    
    # Optional refinements
    cas_variant: CasVariantType = CasVariantType.SPCAS9
    target_region: Optional[str] = None  # "exon", "promoter", "5utr", etc.
    specific_exons: Optional[List[int]] = None
    strand_preference: Optional[Strand] = None
    
    # For base/prime editing
    desired_edit: Optional[str] = None  # e.g., "c.1234A>G" or "p.R123H"
    
    # Context
    cell_type: Optional[str] = None
    delivery_method: Optional[str] = None  # "lentiviral", "RNP", "lipofection"


class ResolvedTarget(BaseModel):
    """Resolved target with full genomic context."""
    specification: TargetSpecification
    chromosome: str
    start: int
    end: int
    strand: Strand
    sequence: str  # Full sequence context
    gene_info: Optional[GeneInfo] = None
    transcript: Optional[Transcript] = None
    annotations: List[Dict[str, Any]] = Field(default_factory=list)


# =============================================================================
# Cas Variants and PAM
# =============================================================================

class PAMVariant(BaseModel):
    """Defines a Cas protein's PAM requirements."""
    name: str
    pam_sequence: str  # IUPAC notation (e.g., "NGG")
    pam_position: PAMPosition
    spacer_length: int
    cut_position: int  # Relative to PAM (negative = upstream)
    
    # Activity modifiers
    pam_scores: Dict[str, float] = Field(default_factory=dict)
    
    # Additional info
    description: Optional[str] = None
    source_organism: Optional[str] = None


# =============================================================================
# Guide Design
# =============================================================================

class SpacerCandidate(BaseModel):
    """A candidate spacer sequence extracted from target region."""
    spacer_sequence: str
    pam_sequence: str
    strand: Strand
    genomic_start: int
    genomic_end: int
    cut_site: int
    context_sequence: str  # Extended context for scoring
    
    # Annotation overlaps
    overlapping_features: List[str] = Field(default_factory=list)
    in_exon: bool = False
    exon_number: Optional[int] = None
    
    @computed_field
    @property
    def full_sequence(self) -> str:
        """Spacer + PAM."""
        return self.spacer_sequence + self.pam_sequence


class EfficiencyScore(BaseModel):
    """Comprehensive efficiency score for a gRNA."""
    overall_score: float = Field(ge=0, le=1)
    
    # Component scores
    components: Dict[str, float] = Field(default_factory=dict)
    
    # Confidence
    confidence_interval: Tuple[float, float] = (0.0, 1.0)
    
    # Interpretation
    interpretation: str = ""
    
    # Model info
    model_version: str = "ensemble_v1"


class DesignedGuide(BaseModel):
    """A fully designed and scored guide RNA."""
    spacer: SpacerCandidate
    efficiency_score: EfficiencyScore
    
    # Composite scores
    composite_score: float = Field(ge=0, le=1)
    position_score: float = Field(ge=0, le=1, default=0.5)
    domain_score: float = Field(ge=0, le=1, default=0.5)
    
    # Rank
    rank: int = 0
    
    # For ordering output
    oligo_forward: Optional[str] = None
    oligo_reverse: Optional[str] = None


# =============================================================================
# Off-Target Analysis
# =============================================================================

class OffTargetSite(BaseModel):
    """An off-target site with risk assessment."""
    chromosome: str
    position: int
    strand: Strand
    sequence: str
    pam: str
    
    # Mismatch info
    mismatches: List[Tuple[int, str, str]] = Field(default_factory=list)  # (pos, ref, alt)
    mismatch_count: int = 0
    
    # Risk scores
    cutting_probability: float = Field(ge=0, le=1)
    risk_score: float = Field(ge=0)
    
    # Annotations
    gene_context: Optional[str] = None
    regulatory_annotation: Optional[str] = None
    in_exon: bool = False
    in_coding_region: bool = False


class RiskDistribution(BaseModel):
    """Probabilistic risk distribution from Monte Carlo simulation."""
    mean: float
    std: float
    percentiles: Dict[int, float] = Field(default_factory=dict)  # {5: x, 25: y, ...}
    n_simulations: int = 10000


class OffTargetAnalysis(BaseModel):
    """Complete off-target analysis results."""
    on_target_spacer: str
    sites: List[OffTargetSite] = Field(default_factory=list)
    
    # Summary stats
    total_sites: int = 0
    sites_in_genes: int = 0
    sites_in_exons: int = 0
    
    # Risk distribution
    risk_distribution: Optional[RiskDistribution] = None
    
    # Analysis metadata
    max_mismatches: int = 4
    genome: Genome = Genome.HUMAN_GRCH38


# =============================================================================
# Red Flags
# =============================================================================

class RedFlag(BaseModel):
    """A warning flag for potential issues."""
    flag_type: str
    severity: RedFlagSeverity
    message: str
    details: Dict[str, Any] = Field(default_factory=dict)


# =============================================================================
# Design Results
# =============================================================================

class DesignResults(BaseModel):
    """Complete results from a guide design run."""
    target: ResolvedTarget
    modality: EditModality
    cas_variant: CasVariantType
    
    # Designed guides
    guides: List[DesignedGuide] = Field(default_factory=list)
    
    # Off-target analysis (optional, computed separately)
    off_target_analyses: Dict[str, OffTargetAnalysis] = Field(default_factory=dict)
    
    # Red flags
    red_flags: List[RedFlag] = Field(default_factory=list)
    
    # Metadata
    created_at: datetime = Field(default_factory=datetime.utcnow)
    scalpel_version: str = "0.1.0"
    parameters: Dict[str, Any] = Field(default_factory=dict)


# =============================================================================
# Experiment Planning
# =============================================================================

class ControlRecommendation(BaseModel):
    """Recommended experimental control."""
    control_type: str  # "positive" or "negative"
    name: str
    description: str
    gene_target: Optional[str] = None
    guide_sequence: Optional[str] = None
    expected_result: str
    interpretation: str


class ValidationAssay(BaseModel):
    """Recommended validation assay."""
    tier: int  # 1, 2, or 3
    name: str
    description: str
    timing: str
    success_criteria: str
    protocol_reference: Optional[str] = None
    estimated_cost: Optional[str] = None
    primers: Optional[Dict[str, str]] = None
    antibodies: Optional[List[str]] = None


class FailureMode(BaseModel):
    """Potential failure mode with troubleshooting."""
    failure: str
    probability: str  # "High", "Medium", "Low"
    symptoms: List[str] = Field(default_factory=list)
    potential_causes: List[str] = Field(default_factory=list)
    troubleshooting: List[str] = Field(default_factory=list)


class ExperimentPlan(BaseModel):
    """Complete experiment planning document."""
    # Overview
    target_gene: str
    edit_modality: EditModality
    
    # Guide information
    selected_guides: List[DesignedGuide] = Field(default_factory=list)
    
    # Controls
    positive_controls: List[ControlRecommendation] = Field(default_factory=list)
    negative_controls: List[ControlRecommendation] = Field(default_factory=list)
    
    # Validation
    validation_assays: List[ValidationAssay] = Field(default_factory=list)
    
    # Failure modes
    failure_modes: List[FailureMode] = Field(default_factory=list)
    
    # Provenance
    created_at: datetime = Field(default_factory=datetime.utcnow)
    scalpel_version: str = "0.1.0"
    input_hash: Optional[str] = None

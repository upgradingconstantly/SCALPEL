"""Models package."""

from scalpel.models.enums import (
    TargetType,
    Genome,
    EditModality,
    CasVariantType,
    PAMPosition,
    Strand,
    RedFlagSeverity,
    ExperimentType,
    ValidationTier,
)
from scalpel.models.data_classes import (
    TargetSpecification,
    PAMVariant,
    SpacerCandidate,
    EfficiencyScore,
    OffTargetSite,
    OffTargetAnalysis,
    RiskDistribution,
    RedFlag,
    DesignedGuide,
    DesignResults,
    GeneInfo,
    Transcript,
    Exon,
    GenomicRegion,
)

__all__ = [
    # Enums
    "TargetType",
    "Genome",
    "EditModality",
    "CasVariantType",
    "PAMPosition",
    "Strand",
    "RedFlagSeverity",
    "ExperimentType",
    "ValidationTier",
    # Data classes
    "TargetSpecification",
    "PAMVariant",
    "SpacerCandidate",
    "EfficiencyScore",
    "OffTargetSite",
    "OffTargetAnalysis",
    "RiskDistribution",
    "RedFlag",
    "DesignedGuide",
    "DesignResults",
    "GeneInfo",
    "Transcript",
    "Exon",
    "GenomicRegion",
]

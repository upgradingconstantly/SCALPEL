"""
SCALPEL - Computational CRISPR Design Platform

A comprehensive platform for gRNA design, off-target analysis,
and experiment planning across multiple CRISPR modalities.
"""

__version__ = "0.1.0"
__author__ = "SCALPEL Team"

from scalpel.models.enums import EditModality, Genome, TargetType
from scalpel.models.data_classes import (
    TargetSpecification,
    SpacerCandidate,
    DesignedGuide,
    EfficiencyScore,
    OffTargetSite,
)

__all__ = [
    # Enums
    "EditModality",
    "Genome", 
    "TargetType",
    # Data classes
    "TargetSpecification",
    "SpacerCandidate",
    "DesignedGuide",
    "EfficiencyScore",
    "OffTargetSite",
]

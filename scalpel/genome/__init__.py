"""Genome data infrastructure package."""

from scalpel.genome.reference import GenomeReference, SequenceUtils
from scalpel.genome.genes import GeneDatabase, get_demo_gene
from scalpel.genome.annotations import AnnotationDatabase, ProteinDomain, RepeatElement
from scalpel.genome.target_resolver import TargetResolver, TargetNotFoundError

__all__ = [
    "GenomeReference",
    "SequenceUtils",
    "GeneDatabase",
    "get_demo_gene",
    "AnnotationDatabase",
    "ProteinDomain",
    "RepeatElement",
    "TargetResolver",
    "TargetNotFoundError",
]

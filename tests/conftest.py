"""
Test configuration and fixtures for SCALPEL.
"""

import pytest
from pathlib import Path
from typing import Generator

from scalpel.models.enums import Genome, EditModality, CasVariantType, Strand
from scalpel.models.data_classes import (
    TargetSpecification,
    TargetType,
    SpacerCandidate,
    GeneInfo,
    Transcript,
    Exon,
)


@pytest.fixture
def sample_sequence() -> str:
    """Sample DNA sequence for testing."""
    return (
        "ATCGATCGATCGATCGATCGATCGGTAGCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGC"
        "GCTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCGCTAGCTAGCTAGCTAGCTAGC"
    )


@pytest.fixture
def sample_spacer_ngg() -> str:
    """Sample 20bp spacer with NGG PAM."""
    return "ATCGATCGATCGATCGATCG"  # Followed by AGG in sequence


@pytest.fixture
def sample_gene_info() -> GeneInfo:
    """Sample gene info for TP53."""
    return GeneInfo(
        gene_id="ENSG00000141510",
        symbol="TP53",
        chromosome="chr17",
        start=7668421,
        end=7687490,
        strand=Strand.MINUS,
        biotype="protein_coding",
        description="tumor protein p53",
        aliases=["p53", "LFS1"],
        transcripts=[
            Transcript(
                transcript_id="ENST00000269305",
                gene_id="ENSG00000141510",
                chromosome="chr17",
                start=7668421,
                end=7687490,
                strand=Strand.MINUS,
                cds_start=7669608,
                cds_end=7676594,
                is_mane_select=True,
                exons=[
                    Exon(exon_number=1, start=7687376, end=7687490, is_constitutive=True),
                    Exon(exon_number=2, start=7676381, end=7676594, is_constitutive=True),
                    Exon(exon_number=3, start=7675993, end=7676272, is_constitutive=True),
                    # ... more exons
                ],
            )
        ],
    )


@pytest.fixture
def sample_target_spec() -> TargetSpecification:
    """Sample target specification."""
    return TargetSpecification(
        target_type=TargetType.GENE_SYMBOL,
        target_value="TP53",
        genome=Genome.HUMAN_GRCH38,
        modality=EditModality.KNOCKOUT,
        cas_variant=CasVariantType.SPCAS9,
    )


@pytest.fixture
def sample_spacer_candidate() -> SpacerCandidate:
    """Sample spacer candidate."""
    return SpacerCandidate(
        spacer_sequence="ATCGATCGATCGATCGATCG",
        pam_sequence="AGG",
        strand=Strand.PLUS,
        genomic_start=100,
        genomic_end=120,
        cut_site=117,
        context_sequence="NNNNNATCGATCGATCGATCGATCGAGGNNNNN",
        in_exon=True,
        exon_number=2,
    )


@pytest.fixture
def temp_data_dir(tmp_path: Path) -> Generator[Path, None, None]:
    """Temporary data directory for tests."""
    data_dir = tmp_path / "scalpel_test"
    data_dir.mkdir()
    (data_dir / "genomes").mkdir()
    (data_dir / "models").mkdir()
    yield data_dir


# Test sequences with known PAM sites
TEST_SEQUENCES = {
    "multiple_ngg": "ATCGATCGATCGATCGATCGAGGCCCGATCGATCGATCGATCGGGGCCCGATCGATCGATCGATCGTGG",
    "no_pam": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    "tttv_pam": "TTTAATCGATCGATCGATCGATCGAATTTCATCGATCGATCGATCGATCGA",  # Cas12a
}

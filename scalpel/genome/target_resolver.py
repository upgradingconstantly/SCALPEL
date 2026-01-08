"""
Target resolution pipeline.

Converts user input (gene symbol, coordinates, sequence) to resolved genomic targets.
"""

from __future__ import annotations

from typing import Optional, List
from dataclasses import dataclass

from scalpel.models.enums import Genome, Strand, TargetType, EditModality
from scalpel.models.data_classes import (
    TargetSpecification,
    ResolvedTarget,
    GeneInfo,
    Transcript,
    GenomicRegion,
)
from scalpel.genome.genes import GeneDatabase, get_demo_gene
from scalpel.genome.reference import GenomeReference


class TargetNotFoundError(Exception):
    """Raised when target cannot be resolved."""
    pass


class TargetResolver:
    """
    Resolves user target specifications to genomic coordinates.
    
    Handles:
    - Gene symbol → coordinates
    - Ensembl ID → coordinates
    - Coordinate parsing
    - Raw sequence input
    - Modality-specific target region selection
    """
    
    def __init__(self, genome: Genome):
        self.genome = genome
        self._gene_db: Optional[GeneDatabase] = None
        self._genome_ref: Optional[GenomeReference] = None
    
    @property
    def gene_db(self) -> GeneDatabase:
        """Lazy-load gene database."""
        if self._gene_db is None:
            self._gene_db = GeneDatabase(self.genome)
        return self._gene_db
    
    @property
    def genome_ref(self) -> Optional[GenomeReference]:
        """Get genome reference if available."""
        # Only return if FASTA file exists
        return self._genome_ref
    
    def resolve(self, spec: TargetSpecification) -> ResolvedTarget:
        """
        Convert any target specification to resolved genomic coordinates.
        
        Args:
            spec: Target specification from user
        
        Returns:
            ResolvedTarget with coordinates, sequence, and annotations
        
        Raises:
            TargetNotFoundError: If target cannot be resolved
        """
        if spec.target_type == TargetType.GENE_SYMBOL:
            return self._resolve_gene_symbol(spec)
        elif spec.target_type == TargetType.ENSEMBL_ID:
            return self._resolve_ensembl(spec)
        elif spec.target_type == TargetType.GENOMIC_COORDINATES:
            return self._resolve_coordinates(spec)
        elif spec.target_type == TargetType.SEQUENCE:
            return self._resolve_sequence(spec)
        else:
            raise ValueError(f"Unknown target type: {spec.target_type}")
    
    def _resolve_gene_symbol(self, spec: TargetSpecification) -> ResolvedTarget:
        """
        Resolve gene symbol to coordinates.
        
        Handles:
        - Symbol aliases
        - Multiple transcripts → canonical selection
        - Modality-specific region selection
        """
        symbol = spec.target_value.upper()
        
        # Try database lookup first
        gene_info = self.gene_db.lookup_symbol(symbol)
        
        # Fall back to demo genes for testing
        if gene_info is None:
            gene_info = get_demo_gene(symbol)
        
        if gene_info is None:
            # Try to suggest alternatives
            suggestions = self.gene_db.suggest_similar(symbol)
            if suggestions:
                raise TargetNotFoundError(
                    f"Gene '{symbol}' not found in {self.genome.value}. "
                    f"Did you mean: {', '.join(suggestions)}?"
                )
            else:
                raise TargetNotFoundError(
                    f"Gene '{symbol}' not found in {self.genome.value}."
                )
        
        # Select transcript
        transcript = self._select_transcript(gene_info, spec)
        
        # Get target region based on modality
        target_region = self._get_target_region(gene_info, transcript, spec)
        
        # Get sequence if genome reference available
        sequence = self._get_sequence(
            gene_info.chromosome,
            target_region.start,
            target_region.end,
            gene_info.strand,
        )
        
        return ResolvedTarget(
            specification=spec,
            chromosome=gene_info.chromosome,
            start=target_region.start,
            end=target_region.end,
            strand=gene_info.strand,
            sequence=sequence,
            gene_info=gene_info,
            transcript=transcript,
        )
    
    def _resolve_ensembl(self, spec: TargetSpecification) -> ResolvedTarget:
        """Resolve Ensembl ID to coordinates."""
        ensembl_id = spec.target_value
        
        gene_info = self.gene_db.lookup_ensembl_id(ensembl_id)
        
        if gene_info is None:
            raise TargetNotFoundError(
                f"Ensembl ID '{ensembl_id}' not found in {self.genome.value}."
            )
        
        transcript = self._select_transcript(gene_info, spec)
        target_region = self._get_target_region(gene_info, transcript, spec)
        sequence = self._get_sequence(
            gene_info.chromosome,
            target_region.start,
            target_region.end,
            gene_info.strand,
        )
        
        return ResolvedTarget(
            specification=spec,
            chromosome=gene_info.chromosome,
            start=target_region.start,
            end=target_region.end,
            strand=gene_info.strand,
            sequence=sequence,
            gene_info=gene_info,
            transcript=transcript,
        )
    
    def _resolve_coordinates(self, spec: TargetSpecification) -> ResolvedTarget:
        """
        Resolve genomic coordinates.
        
        Format: chr17:7668421-7687490 or chr17:7668421-7687490:+
        """
        coord_str = spec.target_value
        
        # Parse format: chr:start-end or chr:start-end:strand
        parts = coord_str.replace(",", "").split(":")
        
        if len(parts) < 2:
            raise ValueError(f"Invalid coordinate format: {coord_str}")
        
        chromosome = parts[0]
        
        if "-" in parts[1]:
            start_str, end_str = parts[1].split("-")
            start = int(start_str)
            end = int(end_str)
        else:
            raise ValueError(f"Invalid coordinate format: {coord_str}")
        
        strand = Strand.PLUS
        if len(parts) >= 3:
            strand = Strand(parts[2])
        
        sequence = self._get_sequence(chromosome, start, end, strand)
        
        return ResolvedTarget(
            specification=spec,
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand,
            sequence=sequence,
        )
    
    def _resolve_sequence(self, spec: TargetSpecification) -> ResolvedTarget:
        """
        Resolve raw sequence input.
        
        Returns target with sequence but no genomic coordinates.
        """
        sequence = spec.target_value.upper()
        
        # Validate sequence
        valid_bases = set("ACGTN")
        if not all(base in valid_bases for base in sequence):
            invalid = set(sequence) - valid_bases
            raise ValueError(f"Invalid bases in sequence: {invalid}")
        
        return ResolvedTarget(
            specification=spec,
            chromosome="unknown",
            start=0,
            end=len(sequence),
            strand=Strand.PLUS,
            sequence=sequence,
        )
    
    def _select_transcript(
        self,
        gene_info: GeneInfo,
        spec: TargetSpecification,
    ) -> Optional[Transcript]:
        """
        Select optimal transcript for targeting.
        
        Priority:
        1. MANE Select transcript
        2. Longest CDS
        3. First available
        """
        transcripts = gene_info.transcripts
        
        if not transcripts:
            return None
        
        # MANE Select is gold standard
        mane = [t for t in transcripts if t.is_mane_select]
        if mane:
            return mane[0]
        
        # Prefer protein coding with longest CDS
        protein_coding = [t for t in transcripts if t.biotype == "protein_coding"]
        if protein_coding:
            return max(protein_coding, key=lambda t: t.cds_length)
        
        return transcripts[0]
    
    def _get_target_region(
        self,
        gene_info: GeneInfo,
        transcript: Optional[Transcript],
        spec: TargetSpecification,
    ) -> GenomicRegion:
        """
        Determine optimal target region based on modality.
        
        Modality-specific logic:
        - KO: Target early exons (avoid last exon → NMD escape)
        - CRISPRi: Target TSS ± 50bp
        - CRISPRa: Target TSS - 400 to -50
        - Base editing: Use full gene region (specific edits handled elsewhere)
        """
        if transcript is None:
            # No transcript info, use full gene region
            return GenomicRegion(
                chromosome=gene_info.chromosome,
                start=gene_info.start,
                end=gene_info.end,
                strand=gene_info.strand,
            )
        
        if spec.modality == EditModality.KNOCKOUT:
            return self._get_knockout_region(gene_info, transcript)
        
        elif spec.modality == EditModality.INTERFERENCE:
            return self._get_crispri_region(gene_info, transcript)
        
        elif spec.modality == EditModality.ACTIVATION:
            return self._get_crispra_region(gene_info, transcript)
        
        else:
            # Base editing, prime editing - use CDS region
            if transcript.cds_start and transcript.cds_end:
                return GenomicRegion(
                    chromosome=gene_info.chromosome,
                    start=transcript.cds_start,
                    end=transcript.cds_end,
                    strand=gene_info.strand,
                )
            return GenomicRegion(
                chromosome=gene_info.chromosome,
                start=gene_info.start,
                end=gene_info.end,
                strand=gene_info.strand,
            )
    
    def _get_knockout_region(
        self,
        gene_info: GeneInfo,
        transcript: Transcript,
    ) -> GenomicRegion:
        """
        Get optimal region for knockout targeting.
        
        Strategy:
        - Skip first exon (alternative start sites)
        - Skip last exon (NMD escape)
        - Prefer early constitutive exons
        """
        exons = transcript.exons
        
        if len(exons) < 3:
            # Short gene - target what we have
            if transcript.cds_start and transcript.cds_end:
                return GenomicRegion(
                    chromosome=gene_info.chromosome,
                    start=transcript.cds_start,
                    end=transcript.cds_end,
                    strand=gene_info.strand,
                )
        
        # Select exons 2-4 (skip first and last)
        target_exons = exons[1:-1] if len(exons) > 2 else exons
        target_exons = target_exons[:3]  # Limit to first 3 valid exons
        
        if not target_exons:
            target_exons = exons
        
        # Get union of selected exons
        start = min(e.start for e in target_exons)
        end = max(e.end for e in target_exons)
        
        return GenomicRegion(
            chromosome=gene_info.chromosome,
            start=start,
            end=end,
            strand=gene_info.strand,
        )
    
    def _get_crispri_region(
        self,
        gene_info: GeneInfo,
        transcript: Transcript,
    ) -> GenomicRegion:
        """
        Get optimal region for CRISPRi targeting.
        
        Target TSS ± 300bp (optimal: TSS to TSS+50).
        """
        tss = transcript.tss
        
        if gene_info.strand == Strand.PLUS:
            start = tss - 50
            end = tss + 300
        else:
            start = tss - 300
            end = tss + 50
        
        return GenomicRegion(
            chromosome=gene_info.chromosome,
            start=max(0, start),
            end=end,
            strand=gene_info.strand,
        )
    
    def _get_crispra_region(
        self,
        gene_info: GeneInfo,
        transcript: Transcript,
    ) -> GenomicRegion:
        """
        Get optimal region for CRISPRa targeting.
        
        Target promoter: TSS -400 to TSS -50.
        """
        tss = transcript.tss
        
        if gene_info.strand == Strand.PLUS:
            start = tss - 400
            end = tss - 50
        else:
            start = tss + 50
            end = tss + 400
        
        return GenomicRegion(
            chromosome=gene_info.chromosome,
            start=max(0, start),
            end=end,
            strand=gene_info.strand,
        )
    
    def _get_sequence(
        self,
        chromosome: str,
        start: int,
        end: int,
        strand: Strand,
    ) -> str:
        """
        Get sequence from genome reference.
        
        Returns placeholder if genome not available.
        """
        if self.genome_ref is not None:
            try:
                return self.genome_ref.get_sequence(chromosome, start, end, strand)
            except Exception:
                pass
        
        # Return placeholder sequence for demo/testing
        # In production, this would require actual genome files
        length = end - start
        return "N" * min(length, 10000)  # Placeholder

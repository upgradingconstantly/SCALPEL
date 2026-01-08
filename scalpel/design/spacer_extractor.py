"""
Spacer extraction engine.

Scans DNA sequences to find PAM sites and extract candidate spacer sequences.
"""

from __future__ import annotations

import re
from typing import List, Optional, Iterator
from dataclasses import dataclass

from scalpel.models.enums import CasVariantType, Strand
from scalpel.models.data_classes import SpacerCandidate
from scalpel.design.cas_variants import (
    PAMVariant,
    CAS_VARIANTS,
    get_cas_variant,
    IUPAC_CODES,
)


class SpacerExtractor:
    """
    Extracts candidate spacer sequences from a target region.
    
    Scans both strands for PAM sites and extracts upstream/downstream
    spacer sequences based on Cas variant requirements.
    """
    
    def __init__(self, cas_variant: CasVariantType = CasVariantType.SPCAS9):
        self.cas_variant = cas_variant
        self.pam_def = get_cas_variant(cas_variant)
    
    def extract_spacers(
        self,
        sequence: str,
        chromosome: str = "unknown",
        start_position: int = 0,
        context_size: int = 30,
    ) -> List[SpacerCandidate]:
        """
        Extract all valid spacer candidates from a sequence.
        
        Args:
            sequence: DNA sequence to scan
            chromosome: Chromosome name for output coordinates
            start_position: Genomic position of sequence start
            context_size: Bases of context to include around spacer
        
        Returns:
            List of SpacerCandidate objects
        """
        sequence = sequence.upper()
        candidates = []
        
        # Scan forward strand
        for candidate in self._scan_strand(
            sequence, chromosome, start_position, Strand.PLUS, context_size
        ):
            candidates.append(candidate)
        
        # Scan reverse strand
        for candidate in self._scan_strand(
            sequence, chromosome, start_position, Strand.MINUS, context_size
        ):
            candidates.append(candidate)
        
        return candidates
    
    def _scan_strand(
        self,
        sequence: str,
        chromosome: str,
        start_position: int,
        strand: Strand,
        context_size: int,
    ) -> Iterator[SpacerCandidate]:
        """
        Scan one strand for PAM sites and extract spacers.
        
        For 3' PAM (Cas9): PAM is after spacer
            Sequence: [spacer][PAM]
        
        For 5' PAM (Cas12a): PAM is before spacer
            Sequence: [PAM][spacer]
        """
        pam_def = self.pam_def
        spacer_len = pam_def.spacer_length
        pam_len = pam_def.pam_length
        
        # Build regex to find PAM sites
        if strand == Strand.PLUS:
            pam_regex = pam_def.pam_regex
        else:
            pam_regex = pam_def.pam_regex_rc
        
        # Determine search sequence
        if strand == Strand.MINUS:
            # For minus strand, we work on reverse complement
            search_seq = self._reverse_complement(sequence)
            # Positions need to be transformed
        else:
            search_seq = sequence
        
        # Find all PAM sites
        for match in pam_regex.finditer(search_seq):
            pam_start = match.start()
            pam_end = match.end()
            pam_seq = match.group()
            
            # Extract spacer based on PAM position
            if pam_def.pam_position.value == "3prime":
                # Spacer is upstream of PAM
                spacer_start = pam_start - spacer_len
                spacer_end = pam_start
            else:
                # Spacer is downstream of PAM (Cas12a)
                spacer_start = pam_end
                spacer_end = pam_end + spacer_len
            
            # Check bounds
            if spacer_start < 0 or spacer_end > len(search_seq):
                continue
            
            spacer_seq = search_seq[spacer_start:spacer_end]
            
            # Validate spacer (no Ns, correct length)
            if not self._is_valid_spacer(spacer_seq):
                continue
            
            # Calculate genomic coordinates
            if strand == Strand.PLUS:
                genomic_start = start_position + spacer_start
                genomic_end = start_position + spacer_end
                genomic_pam_start = start_position + pam_start
            else:
                # Transform coordinates for minus strand
                seq_len = len(sequence)
                genomic_start = start_position + (seq_len - spacer_end)
                genomic_end = start_position + (seq_len - spacer_start)
                genomic_pam_start = start_position + (seq_len - pam_end)
                # Spacer sequence is already reverse complemented from search
            
            # Calculate cut site
            if strand == Strand.PLUS:
                cut_site = start_position + pam_start + pam_def.cut_offset_from_pam
            else:
                cut_site = start_position + (len(sequence) - pam_start) + pam_def.cut_offset_from_pam
            
            # Extract context
            ctx_start = max(0, spacer_start - context_size)
            ctx_end = min(len(search_seq), spacer_end + pam_len + context_size)
            context = search_seq[ctx_start:ctx_end]
            
            yield SpacerCandidate(
                spacer_sequence=spacer_seq,
                pam_sequence=pam_seq,
                strand=strand,
                genomic_start=genomic_start,
                genomic_end=genomic_end,
                cut_site=cut_site,
                context_sequence=context,
            )
    
    def _is_valid_spacer(self, spacer: str) -> bool:
        """Check if spacer is valid (correct length, no Ns)."""
        if len(spacer) != self.pam_def.spacer_length:
            return False
        if "N" in spacer:
            return False
        # Check for valid DNA bases
        if not all(base in "ACGT" for base in spacer):
            return False
        return True
    
    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(complement.get(base, "N") for base in reversed(seq.upper()))


def extract_spacers_from_sequence(
    sequence: str,
    cas_variant: CasVariantType = CasVariantType.SPCAS9,
    chromosome: str = "unknown",
    start_position: int = 0,
) -> List[SpacerCandidate]:
    """
    Convenience function to extract spacers from a sequence.
    
    Args:
        sequence: DNA sequence to scan
        cas_variant: Cas protein variant to use
        chromosome: Chromosome name
        start_position: Genomic start position
    
    Returns:
        List of SpacerCandidate objects
    """
    extractor = SpacerExtractor(cas_variant)
    return extractor.extract_spacers(sequence, chromosome, start_position)


def count_pam_sites(sequence: str, cas_variant: CasVariantType = CasVariantType.SPCAS9) -> int:
    """Quick count of PAM sites in a sequence (both strands)."""
    pam_def = get_cas_variant(cas_variant)
    
    # Forward strand
    forward_count = len(list(pam_def.pam_regex.finditer(sequence.upper())))
    
    # Reverse strand
    rc_seq = SpacerExtractor._reverse_complement(sequence)
    reverse_count = len(list(pam_def.pam_regex.finditer(rc_seq)))
    
    return forward_count + reverse_count

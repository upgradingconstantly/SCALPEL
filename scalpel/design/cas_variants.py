"""
Cas protein variant definitions.

Defines PAM sequences, spacer lengths, and cut positions for all supported Cas variants.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional
from dataclasses import dataclass, field

from scalpel.models.enums import CasVariantType, PAMPosition


# IUPAC nucleotide codes for regex conversion
IUPAC_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",      # Purine
    "Y": "[CT]",      # Pyrimidine
    "S": "[GC]",      # Strong
    "W": "[AT]",      # Weak
    "K": "[GT]",      # Keto
    "M": "[AC]",      # Amino
    "B": "[CGT]",     # Not A
    "D": "[AGT]",     # Not C
    "H": "[ACT]",     # Not G
    "V": "[ACG]",     # Not T
    "N": "[ACGT]",    # Any
}


@dataclass
class PAMVariant:
    """
    Defines a Cas protein's PAM requirements and cleavage properties.
    
    Attributes:
        name: Variant name (e.g., "SpCas9")
        pam_sequence: PAM sequence in IUPAC notation (e.g., "NGG")
        pam_position: Whether PAM is 3' or 5' of spacer
        spacer_length: Length of spacer sequence
        cut_offset_from_pam: Distance from PAM to cut site (negative = upstream)
        staggered_cut: For Cas12a, distance between cuts on opposite strands
        description: Human-readable description
    """
    name: str
    pam_sequence: str
    pam_position: PAMPosition
    spacer_length: int
    cut_offset_from_pam: int
    staggered_cut: int = 0
    description: str = ""
    nickase_available: bool = True
    
    # Cached regex pattern
    _pam_regex: Optional[re.Pattern] = field(default=None, repr=False)
    _pam_regex_rc: Optional[re.Pattern] = field(default=None, repr=False)
    
    @property
    def pam_length(self) -> int:
        """Length of PAM sequence."""
        return len(self.pam_sequence)
    
    @property
    def pam_regex(self) -> re.Pattern:
        """Compiled regex for PAM matching on forward strand."""
        if self._pam_regex is None:
            self._pam_regex = self._compile_pam_regex(self.pam_sequence)
        return self._pam_regex
    
    @property
    def pam_regex_rc(self) -> re.Pattern:
        """Compiled regex for PAM matching on reverse strand."""
        if self._pam_regex_rc is None:
            rc_pam = reverse_complement_iupac(self.pam_sequence)
            self._pam_regex_rc = self._compile_pam_regex(rc_pam)
        return self._pam_regex_rc
    
    def _compile_pam_regex(self, pam: str) -> re.Pattern:
        """Convert IUPAC PAM to regex pattern."""
        regex_parts = []
        for base in pam.upper():
            if base in IUPAC_CODES:
                regex_parts.append(IUPAC_CODES[base])
            else:
                regex_parts.append(base)
        return re.compile("".join(regex_parts))
    
    def matches_pam(self, sequence: str) -> bool:
        """Check if sequence matches PAM pattern."""
        return bool(self.pam_regex.fullmatch(sequence.upper()))
    
    def get_cut_position(self, pam_start: int, strand: str) -> int:
        """
        Calculate cut position given PAM start position.
        
        For 3' PAM (Cas9):
            Spacer is upstream of PAM
            Cut is 3bp upstream of PAM (between positions -4 and -3)
        
        For 5' PAM (Cas12a):
            Spacer is downstream of PAM
            Cut is 18-23bp downstream of PAM end
        
        Args:
            pam_start: 0-based start position of PAM
            strand: "+" or "-"
        
        Returns:
            0-based cut position
        """
        if self.pam_position == PAMPosition.THREE_PRIME:
            # Cas9-style: cut is upstream of PAM
            if strand == "+":
                return pam_start + self.cut_offset_from_pam
            else:
                return pam_start + self.pam_length - self.cut_offset_from_pam
        else:
            # Cas12a-style: cut is downstream of PAM
            if strand == "+":
                return pam_start + self.pam_length + self.cut_offset_from_pam
            else:
                return pam_start - self.cut_offset_from_pam


def reverse_complement_iupac(seq: str) -> str:
    """Reverse complement supporting IUPAC codes."""
    complement = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N",
    }
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))


# =============================================================================
# Pre-defined Cas variants
# =============================================================================

CAS_VARIANTS: Dict[CasVariantType, PAMVariant] = {
    CasVariantType.SPCAS9: PAMVariant(
        name="SpCas9",
        pam_sequence="NGG",
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=20,
        cut_offset_from_pam=-3,  # Cuts 3bp upstream of PAM
        description="S. pyogenes Cas9 - most common, NGG PAM",
    ),
    CasVariantType.SPCAS9_NG: PAMVariant(
        name="SpCas9-NG",
        pam_sequence="NG",
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=20,
        cut_offset_from_pam=-3,
        description="SpCas9 variant with relaxed PAM (NG)",
    ),
    CasVariantType.SPCAS9_VQR: PAMVariant(
        name="SpCas9-VQR",
        pam_sequence="NGA",
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=20,
        cut_offset_from_pam=-3,
        description="SpCas9 variant with NGA PAM",
    ),
    CasVariantType.SACAS9: PAMVariant(
        name="SaCas9",
        pam_sequence="NNGRRT",
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=21,
        cut_offset_from_pam=-3,
        description="S. aureus Cas9 - smaller, NNGRRT PAM",
    ),
    CasVariantType.CAS12A: PAMVariant(
        name="Cas12a",
        pam_sequence="TTTV",
        pam_position=PAMPosition.FIVE_PRIME,
        spacer_length=23,
        cut_offset_from_pam=18,  # Cuts 18bp after PAM
        staggered_cut=5,  # 5bp staggered cut
        description="Cas12a (Cpf1) - 5' PAM, staggered cut",
    ),
    CasVariantType.CAS12A_RR: PAMVariant(
        name="Cas12a-RR",
        pam_sequence="TYCV",
        pam_position=PAMPosition.FIVE_PRIME,
        spacer_length=23,
        cut_offset_from_pam=18,
        staggered_cut=5,
        description="Cas12a RR variant with relaxed PAM",
    ),
    CasVariantType.SPRY: PAMVariant(
        name="SpRY",
        pam_sequence="NRN",  # Also accepts NYN - near PAM-less
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=20,
        cut_offset_from_pam=-3,
        description="SpCas9 variant with near-PAMless activity (NRN/NYN PAM)",
    ),
    CasVariantType.XCAS9: PAMVariant(
        name="xCas9",
        pam_sequence="NG",  # Also accepts GAA, GAT
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=20,
        cut_offset_from_pam=-3,
        description="Expanded PAM recognition (NG, GAA, GAT)",
    ),
    CasVariantType.SACAS9_KKH: PAMVariant(
        name="SaCas9-KKH",
        pam_sequence="NNNRRT",
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=21,
        cut_offset_from_pam=-3,
        description="SaCas9-KKH with relaxed PAM (NNNRRT)",
    ),
    CasVariantType.CJCAS9: PAMVariant(
        name="CjCas9",
        pam_sequence="NNNNRYAC",
        pam_position=PAMPosition.THREE_PRIME,
        spacer_length=22,
        cut_offset_from_pam=-3,
        description="C. jejuni Cas9 - smallest Cas9, good for AAV delivery",
    ),
}


def get_cas_variant(variant_type: CasVariantType) -> PAMVariant:
    """Get PAM variant definition for a Cas type."""
    if variant_type not in CAS_VARIANTS:
        raise ValueError(f"Unknown Cas variant: {variant_type}")
    return CAS_VARIANTS[variant_type]


def list_cas_variants() -> List[CasVariantType]:
    """List all supported Cas variants."""
    return list(CAS_VARIANTS.keys())

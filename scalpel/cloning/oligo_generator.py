"""
Oligo generation for CRISPR cloning workflows.

Generates ready-to-order oligonucleotides for various cloning strategies:
- Golden Gate assembly (BbsI, BsaI, etc.)
- Gibson assembly
- Direct ligation

Output formats:
- IDT ordering format (CSV)
- SnapGene format
- Plain FASTA
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from enum import Enum


class CloningMethod(Enum):
    """Supported cloning methods."""
    GOLDEN_GATE_BBSI = "golden_gate_bbsi"
    GOLDEN_GATE_BSAI = "golden_gate_bsai"
    GIBSON = "gibson"
    DIRECT_LIGATION = "direct_ligation"


@dataclass
class OligoPair:
    """A pair of oligos for cloning a spacer."""
    
    spacer_sequence: str
    forward_name: str
    forward_sequence: str
    reverse_name: str
    reverse_sequence: str
    method: CloningMethod
    notes: str = ""
    
    @property
    def forward_length(self) -> int:
        return len(self.forward_sequence)
    
    @property
    def reverse_length(self) -> int:
        return len(self.reverse_sequence)


# Enzyme-specific overhangs for Golden Gate assembly
GOLDEN_GATE_OVERHANGS = {
    CloningMethod.GOLDEN_GATE_BBSI: {
        "forward_prefix": "CACC",  # BbsI overhang for pSpCas9(BB)-2A-GFP (PX458)
        "reverse_prefix": "AAAC",  # BbsI overhang
        "enzyme": "BbsI",
        "description": "For pSpCas9(BB)-2A-GFP (PX458) and similar vectors",
    },
    CloningMethod.GOLDEN_GATE_BSAI: {
        "forward_prefix": "CACC",  # BsaI overhang
        "reverse_prefix": "AAAC",  # BsaI overhang
        "enzyme": "BsaI",
        "description": "For lentiCRISPRv2 and similar vectors",
    },
}


class OligoGenerator:
    """
    Generate cloning-ready oligonucleotides for CRISPR guides.
    
    Supports multiple cloning strategies and vector systems.
    
    Example:
        >>> generator = OligoGenerator()
        >>> pair = generator.generate_oligos("ATCGATCGATCGATCGATCG", "TP53_g1")
        >>> print(pair.forward_sequence)  # caccATCGATCGATCGATCGATCG
        >>> print(pair.reverse_sequence)  # aaacCGATCGATCGATCGATCGAT
    """
    
    def __init__(
        self,
        method: CloningMethod = CloningMethod.GOLDEN_GATE_BBSI,
        vector_name: Optional[str] = None,
    ):
        """
        Initialize the oligo generator.
        
        Args:
            method: Cloning method to use
            vector_name: Optional vector name for documentation
        """
        self.method = method
        self.vector_name = vector_name or self._default_vector_name()
    
    def _default_vector_name(self) -> str:
        """Get default vector name for the method."""
        defaults = {
            CloningMethod.GOLDEN_GATE_BBSI: "pSpCas9(BB)-2A-GFP (PX458)",
            CloningMethod.GOLDEN_GATE_BSAI: "lentiCRISPRv2",
            CloningMethod.GIBSON: "Custom Gibson vector",
            CloningMethod.DIRECT_LIGATION: "Direct ligation vector",
        }
        return defaults.get(self.method, "Unknown")
    
    def generate_oligos(
        self,
        spacer: str,
        name: str,
        gene: Optional[str] = None,
    ) -> OligoPair:
        """
        Generate a pair of oligos for cloning a spacer.
        
        Args:
            spacer: 20bp spacer sequence (without PAM)
            name: Name prefix for the oligos (e.g., "TP53_g1")
            gene: Optional gene name for documentation
        
        Returns:
            OligoPair with forward and reverse sequences
        """
        spacer = spacer.upper().replace(" ", "")
        
        if self.method in [CloningMethod.GOLDEN_GATE_BBSI, CloningMethod.GOLDEN_GATE_BSAI]:
            return self._generate_golden_gate(spacer, name, gene)
        elif self.method == CloningMethod.GIBSON:
            return self._generate_gibson(spacer, name, gene)
        else:
            return self._generate_direct_ligation(spacer, name, gene)
    
    def _generate_golden_gate(
        self,
        spacer: str,
        name: str,
        gene: Optional[str] = None,
    ) -> OligoPair:
        """Generate Golden Gate assembly oligos."""
        overhangs = GOLDEN_GATE_OVERHANGS[self.method]
        
        # Add G at start if spacer doesn't start with G (for U6 promoter efficiency)
        if not spacer.startswith("G"):
            forward_spacer = "G" + spacer
            needs_g = True
        else:
            forward_spacer = spacer
            needs_g = False
        
        # Forward oligo: overhang + G(optional) + spacer
        forward_seq = overhangs["forward_prefix"].lower() + forward_spacer
        
        # Reverse oligo: overhang + reverse complement of spacer
        reverse_complement = self._reverse_complement(forward_spacer)
        reverse_seq = overhangs["reverse_prefix"].lower() + reverse_complement
        
        notes = f"For {self.vector_name} ({overhangs['enzyme']} digestion)"
        if needs_g:
            notes += ". Added 5' G for U6 promoter efficiency"
        
        return OligoPair(
            spacer_sequence=spacer,
            forward_name=f"{name}_F",
            forward_sequence=forward_seq,
            reverse_name=f"{name}_R",
            reverse_sequence=reverse_seq,
            method=self.method,
            notes=notes,
        )
    
    def _generate_gibson(
        self,
        spacer: str,
        name: str,
        gene: Optional[str] = None,
        homology_length: int = 20,
    ) -> OligoPair:
        """Generate Gibson assembly oligos with homology arms."""
        # Standard U6 promoter 3' end homology
        u6_homology = "GTTTTAGAGCTAGAAATAGC"
        # Standard sgRNA scaffold 5' start homology  
        scaffold_homology = "GTTTAAGAGCTATGCTGGAA"
        
        # Forward: U6 homology + spacer
        forward_seq = u6_homology + spacer
        
        # Reverse: Scaffold homology + reverse complement of spacer
        reverse_complement = self._reverse_complement(spacer)
        reverse_seq = scaffold_homology + reverse_complement
        
        return OligoPair(
            spacer_sequence=spacer,
            forward_name=f"{name}_Gibson_F",
            forward_sequence=forward_seq,
            reverse_name=f"{name}_Gibson_R",
            reverse_sequence=reverse_seq,
            method=self.method,
            notes=f"Gibson assembly oligos with {homology_length}bp homology arms",
        )
    
    def _generate_direct_ligation(
        self,
        spacer: str,
        name: str,
        gene: Optional[str] = None,
    ) -> OligoPair:
        """Generate direct ligation oligos (phosphorylated ends)."""
        forward_seq = spacer
        reverse_seq = self._reverse_complement(spacer)
        
        return OligoPair(
            spacer_sequence=spacer,
            forward_name=f"{name}_F",
            forward_sequence=forward_seq,
            reverse_name=f"{name}_R",
            reverse_sequence=reverse_seq,
            method=self.method,
            notes="Order with 5' phosphorylation for direct ligation",
        )
    
    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(complement.get(base, "N") for base in reversed(sequence.upper()))
    
    def generate_batch(
        self,
        spacers: List[Dict[str, str]],
    ) -> List[OligoPair]:
        """
        Generate oligos for multiple spacers.
        
        Args:
            spacers: List of dicts with 'spacer', 'name', and optional 'gene' keys
        
        Returns:
            List of OligoPair objects
        """
        results = []
        for item in spacers:
            pair = self.generate_oligos(
                spacer=item["spacer"],
                name=item.get("name", f"guide_{len(results)+1}"),
                gene=item.get("gene"),
            )
            results.append(pair)
        return results
    
    def export_idt_format(
        self,
        oligo_pairs: List[OligoPair],
        output_path: Optional[Path] = None,
    ) -> str:
        """
        Export oligos in IDT ordering format (CSV).
        
        Format: Name, Sequence, Scale, Purification
        
        Args:
            oligo_pairs: List of OligoPair objects
            output_path: Optional path to write CSV file
        
        Returns:
            CSV content as string
        """
        output = io.StringIO()
        writer = csv.writer(output)
        
        # IDT format header
        writer.writerow(["Name", "Sequence", "Scale", "Purification"])
        
        for pair in oligo_pairs:
            # Forward oligo
            writer.writerow([
                pair.forward_name,
                pair.forward_sequence,
                "25nm",  # Standard scale
                "STD",   # Standard desalting
            ])
            # Reverse oligo
            writer.writerow([
                pair.reverse_name,
                pair.reverse_sequence,
                "25nm",
                "STD",
            ])
        
        content = output.getvalue()
        
        if output_path:
            output_path.write_text(content)
        
        return content
    
    def export_fasta(
        self,
        oligo_pairs: List[OligoPair],
        output_path: Optional[Path] = None,
    ) -> str:
        """
        Export oligos in FASTA format.
        
        Args:
            oligo_pairs: List of OligoPair objects
            output_path: Optional path to write FASTA file
        
        Returns:
            FASTA content as string
        """
        lines = []
        
        for pair in oligo_pairs:
            lines.append(f">{pair.forward_name}")
            lines.append(pair.forward_sequence)
            lines.append(f">{pair.reverse_name}")
            lines.append(pair.reverse_sequence)
        
        content = "\n".join(lines)
        
        if output_path:
            output_path.write_text(content)
        
        return content


def generate_cloning_oligos(
    spacer: str,
    name: str,
    method: str = "golden_gate_bbsi",
) -> OligoPair:
    """
    Convenience function to generate oligos for a single spacer.
    
    Args:
        spacer: 20bp spacer sequence
        name: Oligo name prefix
        method: Cloning method (golden_gate_bbsi, golden_gate_bsai, gibson, direct_ligation)
    
    Returns:
        OligoPair with forward and reverse sequences
    """
    cloning_method = CloningMethod(method)
    generator = OligoGenerator(method=cloning_method)
    return generator.generate_oligos(spacer, name)

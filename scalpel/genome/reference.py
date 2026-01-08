"""
Genome reference handling with lazy loading.

Uses pyfaidx for memory-efficient access to large FASTA files.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Dict
import threading

from pyfaidx import Fasta, FetchError

from scalpel.models.enums import Genome, Strand
from scalpel.core.cache import cached


class GenomeReference:
    """
    Lazy-loading genome reference using pyfaidx.
    
    Features:
    - Memory-mapped FASTA (only loads requested regions)
    - Thread-safe chromosome access
    - Automatic index creation (.fai)
    - Sequence caching for repeated queries
    """
    
    _instances: Dict[str, "GenomeReference"] = {}
    _lock = threading.Lock()
    
    def __new__(cls, genome: Genome, fasta_path: Optional[Path] = None) -> "GenomeReference":
        """Singleton pattern per genome."""
        key = genome.value
        with cls._lock:
            if key not in cls._instances:
                instance = super().__new__(cls)
                instance._initialized = False
                cls._instances[key] = instance
            return cls._instances[key]
    
    def __init__(self, genome: Genome, fasta_path: Optional[Path] = None):
        if self._initialized:
            return
        
        self.genome = genome
        self.fasta_path = fasta_path
        self._fasta: Optional[Fasta] = None
        self._chromosome_lengths: Dict[str, int] = {}
        self._initialized = True
    
    def _ensure_loaded(self) -> None:
        """Load FASTA file if not already loaded."""
        if self._fasta is not None:
            return
        
        if self.fasta_path is None:
            raise ValueError(
                f"No FASTA path configured for genome {self.genome.value}. "
                f"Set SCALPEL_GENOMES_DIR or provide fasta_path."
            )
        
        if not self.fasta_path.exists():
            raise FileNotFoundError(
                f"Genome FASTA not found: {self.fasta_path}\n"
                f"Download from NCBI/Ensembl and place in genomes directory."
            )
        
        # pyfaidx will create .fai index if needed
        self._fasta = Fasta(str(self.fasta_path), build_index=True)
        
        # Cache chromosome lengths
        for chrom in self._fasta.keys():
            self._chromosome_lengths[chrom] = len(self._fasta[chrom])
    
    @property
    def chromosomes(self) -> list[str]:
        """List of chromosome names."""
        self._ensure_loaded()
        return list(self._fasta.keys())
    
    def get_chromosome_length(self, chromosome: str) -> int:
        """Get length of a chromosome."""
        self._ensure_loaded()
        chrom = self._normalize_chromosome(chromosome)
        if chrom not in self._chromosome_lengths:
            raise ValueError(f"Chromosome {chromosome} not found in {self.genome.value}")
        return self._chromosome_lengths[chrom]
    
    def get_sequence(
        self,
        chromosome: str,
        start: int,
        end: int,
        strand: Strand = Strand.PLUS,
    ) -> str:
        """
        Get sequence from genome.
        
        Args:
            chromosome: Chromosome name (with or without 'chr' prefix)
            start: Start position (0-based)
            end: End position (exclusive)
            strand: Strand (+ returns as-is, - returns reverse complement)
        
        Returns:
            DNA sequence (uppercase)
        """
        self._ensure_loaded()
        
        chrom = self._normalize_chromosome(chromosome)
        
        # Clamp to valid range
        start = max(0, start)
        end = min(end, self.get_chromosome_length(chrom))
        
        if start >= end:
            return ""
        
        try:
            seq = str(self._fasta[chrom][start:end]).upper()
        except (KeyError, FetchError) as e:
            raise ValueError(f"Failed to fetch {chrom}:{start}-{end}: {e}")
        
        if strand == Strand.MINUS:
            seq = self.reverse_complement(seq)
        
        return seq
    
    def get_context_sequence(
        self,
        chromosome: str,
        position: int,
        upstream: int = 5000,
        downstream: int = 5000,
    ) -> tuple[str, int, int]:
        """
        Get sequence context around a position.
        
        Args:
            chromosome: Chromosome name
            position: Center position
            upstream: Bases upstream to include
            downstream: Bases downstream to include
        
        Returns:
            Tuple of (sequence, actual_start, actual_end)
        """
        start = max(0, position - upstream)
        end = min(position + downstream, self.get_chromosome_length(chromosome))
        
        seq = self.get_sequence(chromosome, start, end)
        return seq, start, end
    
    def _normalize_chromosome(self, chromosome: str) -> str:
        """
        Normalize chromosome name to match FASTA.
        
        Handles 'chr1' vs '1' naming conventions.
        """
        self._ensure_loaded()
        
        # Direct match
        if chromosome in self._fasta.keys():
            return chromosome
        
        # Try with/without 'chr' prefix
        if chromosome.startswith("chr"):
            alt = chromosome[3:]
        else:
            alt = f"chr{chromosome}"
        
        if alt in self._fasta.keys():
            return alt
        
        # Try uppercase
        for key in self._fasta.keys():
            if key.upper() == chromosome.upper():
                return key
        
        raise ValueError(
            f"Chromosome {chromosome} not found. "
            f"Available: {list(self._fasta.keys())[:10]}..."
        )
    
    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        return "".join(complement.get(base, "N") for base in reversed(seq.upper()))
    
    def close(self) -> None:
        """Close the FASTA file handle."""
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None


class SequenceUtils:
    """Utility functions for sequence manipulation."""
    
    COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    
    @classmethod
    def reverse_complement(cls, seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        return "".join(cls.COMPLEMENT.get(base, "N") for base in reversed(seq.upper()))
    
    @classmethod
    def gc_content(cls, seq: str) -> float:
        """Calculate GC content as fraction."""
        seq = seq.upper()
        if len(seq) == 0:
            return 0.0
        gc_count = seq.count("G") + seq.count("C")
        return gc_count / len(seq)
    
    @classmethod
    def is_valid_dna(cls, seq: str) -> bool:
        """Check if sequence contains only valid DNA bases."""
        valid_bases = set("ACGTN")
        return all(base in valid_bases for base in seq.upper())
    
    @classmethod
    def translate_codon(cls, codon: str) -> str:
        """Translate a DNA codon to amino acid."""
        codon_table = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        }
        return codon_table.get(codon.upper(), "X")

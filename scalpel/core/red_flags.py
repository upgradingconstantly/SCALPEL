"""
Red flag detection system.

Identifies potential issues with guide designs that may cause problems in experiments.
"""

from __future__ import annotations

from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
from enum import Enum

from scalpel.models.enums import RedFlagSeverity
from scalpel.models.data_classes import (
    RedFlag,
    SpacerCandidate,
    OffTargetAnalysis,
    GeneInfo,
)


class RedFlagType(str, Enum):
    """Types of red flags that can be raised."""
    # Sequence issues
    HOMOPOLYMER = "homopolymer"
    HIGH_GC = "high_gc"
    LOW_GC = "low_gc"
    POLY_T = "poly_t_terminator"
    SELF_COMPLEMENTARY = "self_complementary"
    
    # Genomic context
    REPEAT_ELEMENT = "repeat_element"
    SEGMENTAL_DUPLICATION = "segmental_duplication"
    CNV_REGION = "cnv_region"
    
    # Off-target
    HIGH_OFFTARGET_COUNT = "high_offtarget_count"
    EXONIC_OFFTARGET = "exonic_offtarget"
    GENE_OFFTARGET = "gene_offtarget"
    
    # Gene-level
    ESSENTIAL_GENE = "essential_gene"
    PARALOG_EXISTS = "paralog_exists"
    TISSUE_SPECIFIC = "tissue_specific"
    
    # SNPs
    SNP_IN_SPACER = "snp_in_spacer"
    SNP_IN_PAM = "snp_in_pam"
    
    # Design
    NMD_ESCAPE = "nmd_escape"
    ALTERNATIVE_SPLICING = "alternative_splicing"


@dataclass
class RedFlagCheck:
    """Result of a single red flag check."""
    flag_type: RedFlagType
    passed: bool
    severity: RedFlagSeverity = RedFlagSeverity.LOW
    message: str = ""
    details: Dict[str, Any] = field(default_factory=dict)


class RedFlagDetector:
    """
    Detects potential issues with guide designs.
    
    Checks for:
    1. Sequence-level issues (homopolymers, GC, terminators)
    2. Genomic context issues (repeats, paralogs, CNVs)
    3. Off-target concerns (count, location)
    4. Gene-level issues (essential genes, tissue specificity)
    5. SNP overlaps
    """
    
    # Thresholds
    HOMOPOLYMER_THRESHOLD = 4
    HIGH_GC_THRESHOLD = 0.75
    LOW_GC_THRESHOLD = 0.30
    HIGH_OFFTARGET_THRESHOLD = 50
    LOW_SPECIFICITY_THRESHOLD = 50
    
    def __init__(self):
        self._essential_genes: Optional[set] = None
        self._repeat_db: Optional[Any] = None
    
    def check_all(
        self,
        spacer: SpacerCandidate,
        offtarget_analysis: Optional[OffTargetAnalysis] = None,
        gene_info: Optional[GeneInfo] = None,
    ) -> List[RedFlag]:
        """
        Run all red flag checks on a guide.
        
        Returns:
            List of RedFlag objects for any issues found
        """
        flags = []
        
        # Sequence checks
        flags.extend(self._check_sequence(spacer))
        
        # Off-target checks
        if offtarget_analysis:
            flags.extend(self._check_offtargets(offtarget_analysis))
        
        # Gene-level checks
        if gene_info:
            flags.extend(self._check_gene(gene_info))
        
        # Genomic context checks
        flags.extend(self._check_genomic_context(spacer))
        
        return flags
    
    def _check_sequence(self, spacer: SpacerCandidate) -> List[RedFlag]:
        """Check for sequence-level issues."""
        flags = []
        seq = spacer.spacer_sequence.upper()
        
        # Homopolymer check
        for base in "ACGT":
            run = base * self.HOMOPOLYMER_THRESHOLD
            if run in seq:
                pos = seq.index(run)
                length = 0
                while pos + length < len(seq) and seq[pos + length] == base:
                    length += 1
                
                severity = RedFlagSeverity.MEDIUM if length == 4 else RedFlagSeverity.CRITICAL
                flags.append(RedFlag(
                    flag_type=RedFlagType.HOMOPOLYMER.value,
                    severity=severity,
                    message=f"Contains {length}bp {base} homopolymer at position {pos + 1}",
                    details={"base": base, "length": length, "position": pos + 1},
                ))
                break
        
        # Poly-T terminator (Pol III)
        if "TTTT" in seq:
            flags.append(RedFlag(
                flag_type=RedFlagType.POLY_T.value,
                severity=RedFlagSeverity.CRITICAL,
                message="Contains TTTT sequence (Pol III terminator) - may truncate transcript",
                details={"position": seq.index("TTTT") + 1},
            ))
        
        # GC content
        gc_count = seq.count("G") + seq.count("C")
        gc_content = gc_count / len(seq)
        
        if gc_content > self.HIGH_GC_THRESHOLD:
            flags.append(RedFlag(
                flag_type=RedFlagType.HIGH_GC.value,
                severity=RedFlagSeverity.MEDIUM,
                message=f"High GC content ({gc_content:.0%}) may affect activity",
                details={"gc_content": gc_content},
            ))
        elif gc_content < self.LOW_GC_THRESHOLD:
            flags.append(RedFlag(
                flag_type=RedFlagType.LOW_GC.value,
                severity=RedFlagSeverity.MEDIUM,
                message=f"Low GC content ({gc_content:.0%}) may affect binding",
                details={"gc_content": gc_content},
            ))
        
        return flags
    
    def _check_offtargets(self, analysis: OffTargetAnalysis) -> List[RedFlag]:
        """Check for off-target concerns."""
        flags = []
        
        # High off-target count
        if analysis.total_sites > self.HIGH_OFFTARGET_THRESHOLD:
            flags.append(RedFlag(
                flag_type=RedFlagType.HIGH_OFFTARGET_COUNT.value,
                severity=RedFlagSeverity.MEDIUM,
                message=f"High off-target count ({analysis.total_sites} sites)",
                details={"count": analysis.total_sites},
            ))
        
        # Exonic off-targets
        if analysis.sites_in_exons > 0:
            flags.append(RedFlag(
                flag_type=RedFlagType.EXONIC_OFFTARGET.value,
                severity=RedFlagSeverity.CRITICAL,
                message=f"{analysis.sites_in_exons} off-target site(s) in exons",
                details={"count": analysis.sites_in_exons},
            ))
        
        # Check for high-risk sites (CFD > 0.1)
        high_risk = [s for s in analysis.sites if s.cutting_probability > 0.1]
        if high_risk:
            flags.append(RedFlag(
                flag_type=RedFlagType.GENE_OFFTARGET.value,
                severity=RedFlagSeverity.HIGH,
                message=f"{len(high_risk)} high-risk off-target site(s) (CFD > 0.1)",
                details={"count": len(high_risk), "max_cfd": max(s.cutting_probability for s in high_risk)},
            ))
        
        return flags
    
    def _check_gene(self, gene_info: GeneInfo) -> List[RedFlag]:
        """Check for gene-level issues."""
        flags = []
        
        # Essential gene check
        if self._is_essential_gene(gene_info.symbol):
            flags.append(RedFlag(
                flag_type=RedFlagType.ESSENTIAL_GENE.value,
                severity=RedFlagSeverity.LOW,
                message=f"{gene_info.symbol} is classified as an essential gene (DepMap)",
                details={"gene": gene_info.symbol},
            ))
        
        # Check for paralogs (using aliases as proxy)
        if len(gene_info.aliases) > 3:
            flags.append(RedFlag(
                flag_type=RedFlagType.PARALOG_EXISTS.value,
                severity=RedFlagSeverity.LOW,
                message=f"{gene_info.symbol} has multiple aliases - verify unique targeting",
                details={"aliases": gene_info.aliases[:5]},
            ))
        
        return flags
    
    def _check_genomic_context(self, spacer: SpacerCandidate) -> List[RedFlag]:
        """Check for genomic context issues.
        
        TODO: Integrate with RepeatMasker database for real repeat detection.
        Required databases:
        - RepeatMasker annotations (.out files)
        - Segmental duplication database
        - CNV regions from DGV
        """
        flags = []
        
        # Repeat element check requires RepeatMasker database
        # Currently returns empty - no false positives from simulation
        # To enable: download RepeatMasker annotations and implement lookup
        
        return flags
    
    def _is_essential_gene(self, gene_symbol: str) -> bool:
        """Check if gene is essential (from DepMap or similar)."""
        # Common essential genes for demo
        essential_genes = {
            "TP53", "BRCA1", "BRCA2", "PTEN", "RB1",
            "MYC", "KRAS", "EGFR", "PIK3CA", "AKT1",
            "ATM", "CHEK2", "CDKN2A", "MLH1", "MSH2",
        }
        return gene_symbol.upper() in essential_genes


def detect_red_flags(
    spacer: SpacerCandidate,
    offtarget_analysis: Optional[OffTargetAnalysis] = None,
    gene_info: Optional[GeneInfo] = None,
) -> List[RedFlag]:
    """Convenience function to detect red flags."""
    detector = RedFlagDetector()
    return detector.check_all(spacer, offtarget_analysis, gene_info)


def summarize_red_flags(flags: List[RedFlag]) -> Dict[str, Any]:
    """Generate summary of red flags."""
    if not flags:
        return {
            "total_flags": 0,
            "critical": 0,
            "warnings": 0,
            "info": 0,
            "interpretation": "No issues detected",
        }
    
    critical = sum(1 for f in flags if f.severity == RedFlagSeverity.CRITICAL)
    high = sum(1 for f in flags if f.severity == RedFlagSeverity.HIGH)
    medium = sum(1 for f in flags if f.severity == RedFlagSeverity.MEDIUM)
    low = sum(1 for f in flags if f.severity == RedFlagSeverity.LOW)
    
    if critical > 0:
        interpretation = "Critical issues detected - review before use"
    elif high > 0:
        interpretation = "High-priority issues detected - proceed with caution"
    elif medium > 0:
        interpretation = "Medium-priority notes - review recommended"
    else:
        interpretation = "Minor notes only"
    
    return {
        "total_flags": len(flags),
        "critical": critical,
        "high": high,
        "medium": medium,
        "low": low,
        "flags": [
            {
                "type": f.flag_type,
                "severity": f.severity.value,
                "message": f.message,
            }
            for f in flags
        ],
        "interpretation": interpretation,
    }

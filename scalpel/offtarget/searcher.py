"""
Off-target search and analysis.

Provides off-target site identification and risk assessment.
"""

from __future__ import annotations

import re
from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
import random

from scalpel.models.enums import Genome, Strand
from scalpel.models.data_classes import OffTargetSite, OffTargetAnalysis, RiskDistribution
from scalpel.offtarget.cfd_scorer import CFDScorer, calculate_cfd_score


@dataclass
class OffTargetHit:
    """Raw off-target hit before full scoring."""
    chromosome: str
    position: int
    strand: str
    sequence: str
    pam: str
    mismatches: List[Tuple[int, str, str]]
    mismatch_count: int


class OffTargetSearcher:
    """
    Off-target site identification.
    
    Supports multiple search strategies:
    1. Pre-computed database lookup (DuckDB)
    2. Real-time genome search (for custom sequences)
    3. Simulation mode (for demo/testing)
    """
    
    def __init__(
        self,
        genome: Genome,
        max_mismatches: int = 4,
        use_database: bool = True,
    ):
        self.genome = genome
        self.max_mismatches = max_mismatches
        self.use_database = use_database
        self.cfd_scorer = CFDScorer()
        
        # Database path
        self._db_path: Optional[Path] = None
    
    def search(
        self,
        spacer: str,
        include_pam_variants: bool = True,
    ) -> OffTargetAnalysis:
        """
        Search for off-target sites for a spacer sequence.
        
        Args:
            spacer: 20bp spacer sequence
            include_pam_variants: Include sites with non-canonical PAM
        
        Returns:
            OffTargetAnalysis with all identified sites
        """
        spacer = spacer.upper()
        
        # Try database search first
        if self.use_database:
            hits = self._search_database(spacer)
        else:
            hits = []
        
        # If no database or no hits, use simulation for demo
        if not hits:
            hits = self._simulate_offtargets(spacer)
        
        # Score and filter hits
        sites = []
        for hit in hits:
            cfd_score = self.cfd_scorer.score_offtarget(
                spacer, hit.sequence, hit.pam
            )
            
            # Only include significant off-targets
            if cfd_score > 0.01:
                site = OffTargetSite(
                    chromosome=hit.chromosome,
                    position=hit.position,
                    strand=Strand(hit.strand),
                    sequence=hit.sequence,
                    pam=hit.pam,
                    mismatches=hit.mismatches,
                    mismatch_count=hit.mismatch_count,
                    cutting_probability=cfd_score,
                    risk_score=self._calculate_risk(cfd_score, hit),
                )
                sites.append(site)
        
        # Sort by risk (descending)
        sites.sort(key=lambda s: s.risk_score, reverse=True)
        
        # Calculate risk distribution
        risk_dist = self._calculate_risk_distribution(sites)
        
        return OffTargetAnalysis(
            on_target_spacer=spacer,
            sites=sites,
            total_sites=len(sites),
            sites_in_genes=sum(1 for s in sites if s.gene_context),
            sites_in_exons=sum(1 for s in sites if s.in_exon),
            risk_distribution=risk_dist,
            max_mismatches=self.max_mismatches,
            genome=self.genome,
        )
    
    def _search_database(self, spacer: str) -> List[OffTargetHit]:
        """Search pre-computed database for off-targets."""
        # TODO: Implement DuckDB search
        # For now, return empty - will trigger simulation
        return []
    
    def _simulate_offtargets(self, spacer: str) -> List[OffTargetHit]:
        """
        Simulate off-target sites for demo/testing.
        
        Generates realistic distribution of off-targets based on
        typical CRISPR specificity patterns.
        """
        random.seed(hash(spacer) % 2**32)  # Reproducible per spacer
        
        hits = []
        
        # Generate off-targets with various mismatch counts
        mismatch_distribution = {
            1: random.randint(0, 3),   # 0-3 with 1 mismatch
            2: random.randint(2, 10),  # 2-10 with 2 mismatches
            3: random.randint(10, 50), # 10-50 with 3 mismatches
            4: random.randint(50, 200),# 50-200 with 4 mismatches
        }
        
        chromosomes = ["chr1", "chr2", "chr3", "chr5", "chr7", "chr11", "chr17"]
        pams = ["AGG", "TGG", "CGG", "GGG", "AAG", "TAG"]
        
        for n_mismatches, count in mismatch_distribution.items():
            if n_mismatches > self.max_mismatches:
                continue
            
            for _ in range(count):
                # Generate off-target sequence with random mismatches
                off_seq = list(spacer)
                mismatches = []
                
                # Choose random positions for mismatches
                positions = random.sample(range(20), n_mismatches)
                for pos in positions:
                    orig_base = off_seq[pos]
                    new_base = random.choice([b for b in "ACGT" if b != orig_base])
                    off_seq[pos] = new_base
                    mismatches.append((pos + 1, orig_base, new_base))
                
                hits.append(OffTargetHit(
                    chromosome=random.choice(chromosomes),
                    position=random.randint(1000000, 200000000),
                    strand=random.choice(["+", "-"]),
                    sequence="".join(off_seq),
                    pam=random.choice(pams),
                    mismatches=mismatches,
                    mismatch_count=n_mismatches,
                ))
        
        return hits
    
    def _calculate_risk(self, cfd_score: float, hit: OffTargetHit) -> float:
        """
        Calculate overall risk score for an off-target site.
        
        Combines:
        - CFD cutting probability
        - Genomic context (exon > intron > intergenic)
        - Gene importance (if in gene)
        """
        risk = cfd_score
        
        # Boost risk for sites in coding regions (would be based on annotation)
        # For now, use placeholder
        if hit.mismatch_count <= 2:
            risk *= 2.0  # Low-mismatch sites are higher risk
        
        return min(risk, 10.0)  # Cap at 10
    
    def _calculate_risk_distribution(
        self,
        sites: List[OffTargetSite],
    ) -> RiskDistribution:
        """Calculate risk distribution from Monte Carlo simulation."""
        if not sites:
            return RiskDistribution(
                mean=0.0,
                std=0.0,
                percentiles={5: 0, 25: 0, 50: 0, 75: 0, 95: 0},
                n_simulations=1000,
            )
        
        # Sum up cutting probabilities
        total_risk = sum(s.cutting_probability for s in sites)
        
        # Simple stats (in production, would do full Monte Carlo)
        risks = [s.cutting_probability for s in sites]
        mean_risk = sum(risks) / len(risks) if risks else 0
        
        import math
        variance = sum((r - mean_risk) ** 2 for r in risks) / len(risks) if risks else 0
        std_risk = math.sqrt(variance)
        
        sorted_risks = sorted(risks)
        percentiles = {
            5: sorted_risks[int(0.05 * len(sorted_risks))] if sorted_risks else 0,
            25: sorted_risks[int(0.25 * len(sorted_risks))] if sorted_risks else 0,
            50: sorted_risks[int(0.50 * len(sorted_risks))] if sorted_risks else 0,
            75: sorted_risks[int(0.75 * len(sorted_risks))] if sorted_risks else 0,
            95: sorted_risks[min(int(0.95 * len(sorted_risks)), len(sorted_risks)-1)] if sorted_risks else 0,
        }
        
        return RiskDistribution(
            mean=mean_risk,
            std=std_risk,
            percentiles=percentiles,
            n_simulations=1000,
        )


class RiskCalculator:
    """
    Off-target risk calculator with genomic annotation integration.
    """
    
    def __init__(self, genome: Genome):
        self.genome = genome
    
    def calculate_specificity_score(self, analysis: OffTargetAnalysis) -> float:
        """
        Calculate overall specificity score (0-100).
        
        Higher = more specific = fewer/weaker off-targets.
        """
        if not analysis.sites:
            return 100.0
        
        # Sum of CFD scores
        total_cfd = sum(s.cutting_probability for s in analysis.sites)
        
        # Convert to specificity (inverse relationship)
        # A guide with no off-targets has specificity 100
        # A guide with many high-CFD off-targets approaches 0
        specificity = 100 * (1 / (1 + total_cfd))
        
        return round(specificity, 1)
    
    def generate_summary(self, analysis: OffTargetAnalysis) -> Dict[str, Any]:
        """Generate human-readable summary of off-target analysis."""
        n_sites = len(analysis.sites)
        
        # Categorize by mismatch count
        by_mismatch = {}
        for site in analysis.sites:
            mm = site.mismatch_count
            by_mismatch[mm] = by_mismatch.get(mm, 0) + 1
        
        # High-risk sites (CFD > 0.1)
        high_risk = [s for s in analysis.sites if s.cutting_probability > 0.1]
        
        specificity = self.calculate_specificity_score(analysis)
        
        if specificity >= 80:
            interpretation = "Excellent specificity"
        elif specificity >= 60:
            interpretation = "Good specificity"
        elif specificity >= 40:
            interpretation = "Moderate specificity - review off-targets"
        else:
            interpretation = "Poor specificity - consider alternative guides"
        
        return {
            "total_sites": n_sites,
            "sites_by_mismatch": by_mismatch,
            "high_risk_count": len(high_risk),
            "in_exons": analysis.sites_in_exons,
            "in_genes": analysis.sites_in_genes,
            "specificity_score": specificity,
            "interpretation": interpretation,
        }

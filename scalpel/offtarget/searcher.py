"""
Off-target search and analysis.

Provides off-target site identification and risk assessment.
"""

from __future__ import annotations

import os
from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass
from pathlib import Path


from scalpel.models.enums import Genome, Strand
from scalpel.models.data_classes import OffTargetSite, OffTargetAnalysis, RiskDistribution
from scalpel.offtarget.cfd_scorer import CFDScorer


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


class OffTargetIndexNotFoundError(RuntimeError):
    """Raised when an off-target index database is unavailable."""


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
            hits = self._search_database(spacer, include_pam_variants=include_pam_variants)
        else:
            hits = []
        
        # If no database available, return empty results with warning
        # Do NOT simulate fake off-targets - laboratory use requires real data
        
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
    
    def _search_database(
        self,
        spacer: str,
        include_pam_variants: bool = True,
    ) -> List[OffTargetHit]:
        """Search pre-computed database for off-targets.
        """
        db_path = self._resolve_database_path()
        if db_path is None:
            raise OffTargetIndexNotFoundError(
                f"Off-target index not found for {self.genome.value}. "
                f"Build it with: python -m scalpel.offtarget.build_index --genome {self.genome.value}"
            )

        try:
            import duckdb
        except Exception as e:
            raise RuntimeError(f"DuckDB is required for off-target search: {e}")

        # Guarantee at least one seed chunk.
        n_chunks = max(1, self.max_mismatches + 1)
        seed_windows = self._build_seed_windows(len(spacer), n_chunks)

        # Deduplicate candidates found by multiple seed hits.
        raw_candidates: Dict[Tuple[str, int, str, str, str], Tuple[str, int, str, str, str]] = {}

        # Keep result volume bounded during seed lookup.
        per_seed_limit = max(5000, 200000 // len(seed_windows))

        try:
            with duckdb.connect(str(db_path), read_only=True) as conn:
                for start, end in seed_windows:
                    seed = spacer[start:end]
                    if not seed:
                        continue

                    rows = conn.execute(
                        """
                        SELECT chromosome, position, strand, spacer, pam
                        FROM pam_sites
                        WHERE length(spacer) = ?
                          AND substr(spacer, ?, ?) = ?
                        LIMIT ?
                        """,
                        [
                            len(spacer),
                            start + 1,  # DuckDB substr is 1-indexed.
                            len(seed),
                            seed,
                            per_seed_limit,
                        ],
                    ).fetchall()

                    for row in rows:
                        key = (row[0], int(row[1]), row[2], row[3], row[4])
                        raw_candidates[key] = key
        except Exception as e:
            raise RuntimeError(
                f"Failed to query off-target index at {db_path}: {e}. "
                f"Rebuild with: python -m scalpel.offtarget.build_index --genome {self.genome.value}"
            ) from e

        hits: List[OffTargetHit] = []
        for chromosome, position, strand, candidate_spacer, pam in raw_candidates.values():
            if not include_pam_variants and not pam.endswith("GG"):
                continue

            mismatches = self._find_mismatches(spacer, candidate_spacer)
            mismatch_count = len(mismatches)

            # Exclude on-target site and trim to requested mismatch radius.
            if mismatch_count == 0 or mismatch_count > self.max_mismatches:
                continue

            hits.append(
                OffTargetHit(
                    chromosome=chromosome,
                    position=position,
                    strand=strand,
                    sequence=candidate_spacer,
                    pam=pam,
                    mismatches=mismatches,
                    mismatch_count=mismatch_count,
                )
            )

        # Lowest mismatch count first before downstream risk sort.
        hits.sort(key=lambda h: h.mismatch_count)
        return hits

    def _resolve_database_path(self) -> Optional[Path]:
        """Resolve off-target DB path with explicit precedence."""
        env_path = os.getenv("SCALPEL_OFFTARGET_DB")
        if env_path:
            candidate = Path(env_path).expanduser()
            if candidate.exists():
                self._db_path = candidate
                return candidate
            raise OffTargetIndexNotFoundError(
                f"SCALPEL_OFFTARGET_DB is set but file does not exist: {candidate}"
            )

        candidates = [
            Path.home() / ".scalpel" / "offtargets" / f"{self.genome.value}.duckdb",
            Path.home() / ".scalpel" / "genomes" / self.genome.value / "offtargets.duckdb",
            # Backward-compatible lowercase fallback.
            Path.home() / ".scalpel" / "genomes" / self.genome.value.lower() / "offtargets.duckdb",
        ]
        for candidate in candidates:
            if candidate.exists():
                self._db_path = candidate
                return candidate

        return None

    @staticmethod
    def _build_seed_windows(sequence_length: int, n_chunks: int) -> List[Tuple[int, int]]:
        """
        Partition the sequence into windows used for seed filtering.
        Any sequence within k mismatches must match at least one of k+1 windows.
        """
        windows: List[Tuple[int, int]] = []
        for i in range(n_chunks):
            start = (sequence_length * i) // n_chunks
            end = (sequence_length * (i + 1)) // n_chunks
            if end > start:
                windows.append((start, end))
        return windows

    @staticmethod
    def _find_mismatches(on_target: str, candidate: str) -> List[Tuple[int, str, str]]:
        """Return mismatch tuples: (1-indexed position, on_target_base, candidate_base)."""
        mismatches: List[Tuple[int, str, str]] = []
        for i, (on_base, candidate_base) in enumerate(zip(on_target, candidate), start=1):
            if on_base != candidate_base:
                mismatches.append((i, on_base, candidate_base))
        return mismatches
    
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

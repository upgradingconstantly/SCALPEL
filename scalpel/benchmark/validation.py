"""
Benchmark and validation suite for SCALPEL.

Provides tools for evaluating gRNA design performance against published datasets.
"""

from __future__ import annotations

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field
import math


@dataclass
class BenchmarkResult:
    """Result from a single benchmark run."""
    benchmark_name: str
    n_guides: int
    
    # Correlation metrics
    spearman_correlation: float
    pearson_correlation: float
    
    # Classification metrics (if binary outcomes available)
    auroc: Optional[float] = None
    auprc: Optional[float] = None
    
    # Additional stats
    mean_predicted: float = 0.0
    mean_actual: float = 0.0
    rmse: float = 0.0


@dataclass 
class BenchmarkDataset:
    """A benchmark dataset for validation."""
    name: str
    description: str
    source: str
    n_guides: int
    spacers: List[str] = field(default_factory=list)
    efficiency_scores: List[float] = field(default_factory=list)
    has_binary_outcome: bool = False


class BenchmarkSuite:
    """
    Benchmark suite for evaluating gRNA efficiency predictions.
    
    Supports validation against published datasets:
    - Doench et al. 2016 (Rule Set 2)
    - Wang et al. 2019
    - Kim et al. 2020
    """
    
    # Known benchmark datasets
    DATASETS = {
        "doench2016": BenchmarkDataset(
            name="Doench 2016",
            description="Rule Set 2 training data",
            source="doi:10.1038/nbt.3437",
            n_guides=2076,
        ),
        "wang2019": BenchmarkDataset(
            name="Wang 2019", 
            description="Deep learning training set",
            source="doi:10.1038/s41587-019-0095-1",
            n_guides=15000,
        ),
    }
    
    def __init__(self):
        self.results: List[BenchmarkResult] = []
    
    def run_benchmark(
        self,
        dataset_name: str,
        predictions: List[float],
        actuals: List[float],
    ) -> BenchmarkResult:
        """
        Run benchmark comparison between predictions and actuals.
        
        Args:
            dataset_name: Name of the benchmark dataset
            predictions: Predicted efficiency scores
            actuals: Actual measured efficiencies
            
        Returns:
            BenchmarkResult with correlation metrics
        """
        if len(predictions) != len(actuals):
            raise ValueError("Predictions and actuals must have same length")
        
        n = len(predictions)
        
        # Calculate Spearman correlation
        spearman = self._spearman_correlation(predictions, actuals)
        
        # Calculate Pearson correlation
        pearson = self._pearson_correlation(predictions, actuals)
        
        # Calculate RMSE
        rmse = math.sqrt(sum((p - a) ** 2 for p, a in zip(predictions, actuals)) / n)
        
        result = BenchmarkResult(
            benchmark_name=dataset_name,
            n_guides=n,
            spearman_correlation=spearman,
            pearson_correlation=pearson,
            mean_predicted=sum(predictions) / n,
            mean_actual=sum(actuals) / n,
            rmse=rmse,
        )
        
        self.results.append(result)
        return result
    
    def _spearman_correlation(self, x: List[float], y: List[float]) -> float:
        """Calculate Spearman rank correlation."""
        n = len(x)
        if n < 2:
            return 0.0
        
        # Compute ranks
        x_ranks = self._compute_ranks(x)
        y_ranks = self._compute_ranks(y)
        
        # Spearman = Pearson of ranks
        return self._pearson_correlation(x_ranks, y_ranks)
    
    def _pearson_correlation(self, x: List[float], y: List[float]) -> float:
        """Calculate Pearson correlation coefficient."""
        n = len(x)
        if n < 2:
            return 0.0
        
        mean_x = sum(x) / n
        mean_y = sum(y) / n
        
        numerator = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
        
        sum_sq_x = sum((xi - mean_x) ** 2 for xi in x)
        sum_sq_y = sum((yi - mean_y) ** 2 for yi in y)
        
        denominator = math.sqrt(sum_sq_x * sum_sq_y)
        
        if denominator == 0:
            return 0.0
        
        return numerator / denominator
    
    def _compute_ranks(self, values: List[float]) -> List[float]:
        """Compute ranks for a list of values (handles ties with average rank)."""
        n = len(values)
        indexed = [(v, i) for i, v in enumerate(values)]
        sorted_indexed = sorted(indexed, key=lambda x: x[0])
        
        ranks = [0.0] * n
        i = 0
        while i < n:
            # Find ties
            j = i
            while j < n and sorted_indexed[j][0] == sorted_indexed[i][0]:
                j += 1
            
            # Average rank for ties
            avg_rank = (i + j + 1) / 2
            for k in range(i, j):
                ranks[sorted_indexed[k][1]] = avg_rank
            
            i = j
        
        return ranks
    
    def generate_report(self) -> Dict[str, Any]:
        """Generate summary report of all benchmark results."""
        if not self.results:
            return {"message": "No benchmarks run yet"}
        
        return {
            "n_benchmarks": len(self.results),
            "results": [
                {
                    "dataset": r.benchmark_name,
                    "n_guides": r.n_guides,
                    "spearman": round(r.spearman_correlation, 3),
                    "pearson": round(r.pearson_correlation, 3),
                    "rmse": round(r.rmse, 4),
                }
                for r in self.results
            ],
            "average_spearman": round(
                sum(r.spearman_correlation for r in self.results) / len(self.results), 3
            ),
        }


def run_quick_validation() -> BenchmarkResult:
    """
    Run a quick validation with synthetic data.
    
    Useful for testing the benchmark pipeline.
    """
    suite = BenchmarkSuite()
    
    # Synthetic test data (simulates decent correlation)
    predictions = [0.3, 0.5, 0.7, 0.4, 0.8, 0.6, 0.2, 0.9, 0.5, 0.7]
    actuals = [0.35, 0.45, 0.75, 0.42, 0.82, 0.55, 0.28, 0.85, 0.52, 0.68]
    
    return suite.run_benchmark("synthetic_test", predictions, actuals)

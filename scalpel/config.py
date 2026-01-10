"""
Configuration system for SCALPEL.

Uses Pydantic Settings for type-safe configuration with environment variable support.
"""

from pathlib import Path
from typing import Optional, Dict, Any
from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict
from functools import lru_cache


class GenomeConfig(BaseSettings):
    """Configuration for a specific genome."""
    name: str
    fasta_path: Optional[Path] = None
    fai_path: Optional[Path] = None
    annotation_path: Optional[Path] = None
    precomputed_ot_path: Optional[Path] = None


class CacheConfig(BaseSettings):
    """Cache configuration."""
    # L1: In-memory LRU
    memory_cache_size: int = 10000
    
    # L2: Redis (optional)
    redis_enabled: bool = False
    redis_url: str = "redis://localhost:6379/0"
    redis_ttl: int = 3600  # seconds
    
    # L3: SQLite/DuckDB
    db_cache_path: Path = Path("~/.scalpel/cache.db").expanduser()


class ModelConfig(BaseSettings):
    """ML model configuration."""
    efficiency_model_path: Optional[Path] = None
    use_onnx: bool = True
    device: str = "cpu"  # "cpu", "cuda", "mps"
    batch_size: int = 32


class DesignConfig(BaseSettings):
    """Design algorithm configuration."""
    # Spacer extraction
    default_cas_variant: str = "SpCas9"
    
    # Efficiency scoring weights
    efficiency_weights: Dict[str, float] = Field(default_factory=lambda: {
        "gc_content": 0.10,
        "gc_in_seed": 0.10,
        "homopolymer": 0.10,
        "position_bias": 0.05,
        "secondary_structure": 0.10,
        "rule_set_2": 0.25,
        "ml_score": 0.30,
    })
    
    # Off-target search
    max_mismatches: int = 4
    include_pam_mismatches: bool = True
    
    # Output
    default_n_guides: int = 10


class SCALPELConfig(BaseSettings):
    """Main configuration for SCALPEL."""
    
    model_config = SettingsConfigDict(
        env_prefix="SCALPEL_",
        env_nested_delimiter="__",
        extra="ignore",
    )
    
    # Paths
    data_dir: Path = Path("~/.scalpel").expanduser()
    genomes_dir: Path = Path("~/.scalpel/genomes").expanduser()
    models_dir: Path = Path("~/.scalpel/models").expanduser()
    
    # Sub-configs
    cache: CacheConfig = Field(default_factory=CacheConfig)
    models: ModelConfig = Field(default_factory=ModelConfig)
    design: DesignConfig = Field(default_factory=DesignConfig)
    
    # Logging
    log_level: str = "INFO"
    log_file: Optional[Path] = None
    
    # API
    api_host: str = "127.0.0.1"
    api_port: int = 8000
    
    def ensure_directories(self) -> None:
        """Create required directories if they don't exist."""
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.models_dir.mkdir(parents=True, exist_ok=True)
    
    def get_genome_config(self, genome_name: str) -> GenomeConfig:
        """Get configuration for a specific genome.
        
        Searches multiple locations:
        1. ~/.scalpel/genomes/{genome}/
        2. Project-local scalpel/genomes/
        """
        # Try standard location first
        genome_dir = self.genomes_dir / genome_name.lower()
        
        # Also check project-local genomes directory
        project_genomes_dir = Path(__file__).parent / "genomes"
        
        # Find FASTA file - check multiple naming conventions
        fasta_path = None
        possible_names = [
            f"{genome_name}.fa",
            f"{genome_name}.fasta",
            # Human genomes (Ensembl naming)
            f"Homo_sapiens.{genome_name}.dna.primary_assembly.fa",
            f"Homo_sapiens.{genome_name}.dna.primary_assembly.fasta",
            # Mouse genomes (Ensembl naming)
            f"Mus_musculus.{genome_name}.dna.primary_assembly.fa",
            f"Mus_musculus.{genome_name}.dna.primary_assembly.fasta",
        ]
        
        # Check standard location
        if genome_dir.exists():
            for name in possible_names:
                candidate = genome_dir / name
                if candidate.exists():
                    fasta_path = candidate
                    break
        
        # Check project-local location if not found
        if fasta_path is None and project_genomes_dir.exists():
            for name in possible_names:
                candidate = project_genomes_dir / name
                if candidate.exists():
                    fasta_path = candidate
                    break
        
        return GenomeConfig(
            name=genome_name,
            fasta_path=fasta_path,
            fai_path=fasta_path.with_suffix(".fa.fai") if fasta_path else None,
            annotation_path=genome_dir / "annotations.gtf" if genome_dir.exists() else None,
            precomputed_ot_path=genome_dir / "offtargets.duckdb" if genome_dir.exists() else None,
        )


@lru_cache()
def get_config() -> SCALPELConfig:
    """Get cached configuration singleton."""
    config = SCALPELConfig()
    config.ensure_directories()
    return config


# Convenience accessor
config = get_config()

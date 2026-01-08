"""
Gene database interface.

Provides gene lookup, transcript selection, and annotation caching.
Uses SQLite for local caching of gene annotations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, List, Dict, Any
import sqlite3
import json
from dataclasses import dataclass

from scalpel.models.enums import Genome, Strand
from scalpel.models.data_classes import GeneInfo, Transcript, Exon
from scalpel.core.cache import cached


class GeneDatabase:
    """
    Gene annotation database with caching.
    
    Features:
    - SQLite-backed local cache
    - Gene symbol and alias resolution
    - Transcript selection (MANE Select priority)
    - Exon boundary retrieval
    """
    
    def __init__(self, genome: Genome, db_path: Optional[Path] = None):
        self.genome = genome
        
        if db_path is None:
            from scalpel.config import get_config
            config = get_config()
            db_path = config.data_dir / "genes" / f"{genome.value}.db"
        
        self.db_path = db_path
        self._connection: Optional[sqlite3.Connection] = None
        self._ensure_db()
    
    def _ensure_db(self) -> None:
        """Create database and tables if needed."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        conn = self._get_connection()
        conn.executescript("""
            CREATE TABLE IF NOT EXISTS genes (
                gene_id TEXT PRIMARY KEY,
                symbol TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                strand TEXT NOT NULL,
                biotype TEXT,
                description TEXT,
                data JSON
            );
            
            CREATE TABLE IF NOT EXISTS aliases (
                alias TEXT PRIMARY KEY,
                gene_id TEXT NOT NULL,
                FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
            );
            
            CREATE TABLE IF NOT EXISTS transcripts (
                transcript_id TEXT PRIMARY KEY,
                gene_id TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                strand TEXT NOT NULL,
                cds_start INTEGER,
                cds_end INTEGER,
                biotype TEXT,
                is_mane_select INTEGER DEFAULT 0,
                data JSON,
                FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
            );
            
            CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol);
            CREATE INDEX IF NOT EXISTS idx_genes_chromosome ON genes(chromosome, start, end);
            CREATE INDEX IF NOT EXISTS idx_transcripts_gene ON transcripts(gene_id);
        """)
        conn.commit()
    
    def _get_connection(self) -> sqlite3.Connection:
        """Get database connection."""
        if self._connection is None:
            self._connection = sqlite3.connect(str(self.db_path))
            self._connection.row_factory = sqlite3.Row
        return self._connection
    
    def lookup_symbol(self, symbol: str) -> Optional[GeneInfo]:
        """
        Lookup gene by symbol.
        
        Args:
            symbol: Gene symbol (e.g., "TP53", "BRCA1")
        
        Returns:
            GeneInfo if found, None otherwise
        """
        conn = self._get_connection()
        
        # Try exact match first (case-insensitive)
        cursor = conn.execute(
            "SELECT * FROM genes WHERE UPPER(symbol) = UPPER(?)",
            (symbol,)
        )
        row = cursor.fetchone()
        
        if row:
            return self._row_to_gene_info(row)
        
        # Try alias lookup
        return self.lookup_alias(symbol)
    
    def lookup_alias(self, alias: str) -> Optional[GeneInfo]:
        """Lookup gene by alias."""
        conn = self._get_connection()
        
        cursor = conn.execute(
            """
            SELECT g.* FROM genes g
            JOIN aliases a ON g.gene_id = a.gene_id
            WHERE UPPER(a.alias) = UPPER(?)
            """,
            (alias,)
        )
        row = cursor.fetchone()
        
        if row:
            return self._row_to_gene_info(row)
        return None
    
    def lookup_ensembl_id(self, ensembl_id: str) -> Optional[GeneInfo]:
        """Lookup gene by Ensembl ID."""
        conn = self._get_connection()
        
        # Handle version suffix (ENSG00000141510.16 -> ENSG00000141510)
        base_id = ensembl_id.split(".")[0]
        
        cursor = conn.execute(
            "SELECT * FROM genes WHERE gene_id = ? OR gene_id LIKE ?",
            (base_id, f"{base_id}.%")
        )
        row = cursor.fetchone()
        
        if row:
            return self._row_to_gene_info(row)
        return None
    
    def get_genes_at_position(
        self,
        chromosome: str,
        start: int,
        end: int,
    ) -> List[GeneInfo]:
        """Get genes overlapping a genomic region."""
        conn = self._get_connection()
        
        cursor = conn.execute(
            """
            SELECT * FROM genes
            WHERE chromosome = ? AND start < ? AND end > ?
            ORDER BY start
            """,
            (chromosome, end, start)
        )
        
        return [self._row_to_gene_info(row) for row in cursor.fetchall()]
    
    def suggest_similar(self, symbol: str, limit: int = 5) -> List[str]:
        """Suggest similar gene symbols for typo correction."""
        conn = self._get_connection()
        
        # Simple prefix matching
        cursor = conn.execute(
            """
            SELECT DISTINCT symbol FROM genes
            WHERE symbol LIKE ? OR symbol LIKE ?
            LIMIT ?
            """,
            (f"{symbol}%", f"%{symbol}%", limit)
        )
        
        return [row["symbol"] for row in cursor.fetchall()]
    
    def get_transcripts(self, gene_id: str) -> List[Transcript]:
        """Get all transcripts for a gene."""
        conn = self._get_connection()
        
        cursor = conn.execute(
            "SELECT * FROM transcripts WHERE gene_id = ? ORDER BY is_mane_select DESC, end - start DESC",
            (gene_id,)
        )
        
        transcripts = []
        for row in cursor.fetchall():
            data = json.loads(row["data"]) if row["data"] else {}
            exons = [Exon(**e) for e in data.get("exons", [])]
            
            transcripts.append(Transcript(
                transcript_id=row["transcript_id"],
                gene_id=row["gene_id"],
                chromosome=row["chromosome"],
                start=row["start"],
                end=row["end"],
                strand=Strand(row["strand"]),
                cds_start=row["cds_start"],
                cds_end=row["cds_end"],
                biotype=row["biotype"] or "protein_coding",
                is_mane_select=bool(row["is_mane_select"]),
                exons=exons,
            ))
        
        return transcripts
    
    def _row_to_gene_info(self, row: sqlite3.Row) -> GeneInfo:
        """Convert database row to GeneInfo object."""
        data = json.loads(row["data"]) if row["data"] else {}
        
        gene_info = GeneInfo(
            gene_id=row["gene_id"],
            symbol=row["symbol"],
            chromosome=row["chromosome"],
            start=row["start"],
            end=row["end"],
            strand=Strand(row["strand"]),
            biotype=row["biotype"] or "protein_coding",
            description=row["description"],
            aliases=data.get("aliases", []),
        )
        
        # Load transcripts
        gene_info.transcripts = self.get_transcripts(row["gene_id"])
        
        return gene_info
    
    # ==========================================================================
    # Data import methods
    # ==========================================================================
    
    def import_gene(self, gene: GeneInfo) -> None:
        """Import a gene into the database."""
        conn = self._get_connection()
        
        data = {
            "aliases": gene.aliases,
        }
        
        conn.execute(
            """
            INSERT OR REPLACE INTO genes 
            (gene_id, symbol, chromosome, start, end, strand, biotype, description, data)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                gene.gene_id,
                gene.symbol,
                gene.chromosome,
                gene.start,
                gene.end,
                gene.strand.value,
                gene.biotype,
                gene.description,
                json.dumps(data),
            )
        )
        
        # Import aliases
        for alias in gene.aliases:
            conn.execute(
                "INSERT OR IGNORE INTO aliases (alias, gene_id) VALUES (?, ?)",
                (alias, gene.gene_id)
            )
        
        # Import transcripts
        for transcript in gene.transcripts:
            self.import_transcript(transcript)
        
        conn.commit()
    
    def import_transcript(self, transcript: Transcript) -> None:
        """Import a transcript into the database."""
        conn = self._get_connection()
        
        data = {
            "exons": [
                {
                    "exon_number": e.exon_number,
                    "start": e.start,
                    "end": e.end,
                    "is_constitutive": e.is_constitutive,
                }
                for e in transcript.exons
            ]
        }
        
        conn.execute(
            """
            INSERT OR REPLACE INTO transcripts
            (transcript_id, gene_id, chromosome, start, end, strand, 
             cds_start, cds_end, biotype, is_mane_select, data)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                transcript.transcript_id,
                transcript.gene_id,
                transcript.chromosome,
                transcript.start,
                transcript.end,
                transcript.strand.value,
                transcript.cds_start,
                transcript.cds_end,
                transcript.biotype,
                1 if transcript.is_mane_select else 0,
                json.dumps(data),
            )
        )
    
    def close(self) -> None:
        """Close database connection."""
        if self._connection:
            self._connection.close()
            self._connection = None


# =============================================================================
# Built-in gene data for demo/testing
# =============================================================================

DEMO_GENES: Dict[str, Dict[str, Any]] = {
    "TP53": {
        "gene_id": "ENSG00000141510",
        "symbol": "TP53",
        "chromosome": "chr17",
        "start": 7668421,
        "end": 7687490,
        "strand": "-",
        "biotype": "protein_coding",
        "description": "tumor protein p53",
        "aliases": ["p53", "LFS1", "TRP53"],
        "transcripts": [
            {
                "transcript_id": "ENST00000269305",
                "cds_start": 7669608,
                "cds_end": 7676594,
                "is_mane_select": True,
                "exons": [
                    {"exon_number": 1, "start": 7687376, "end": 7687490},
                    {"exon_number": 2, "start": 7676381, "end": 7676594},
                    {"exon_number": 3, "start": 7675993, "end": 7676272},
                    {"exon_number": 4, "start": 7675052, "end": 7675236},
                    {"exon_number": 5, "start": 7674858, "end": 7674970},
                    {"exon_number": 6, "start": 7674180, "end": 7674289},
                    {"exon_number": 7, "start": 7673700, "end": 7673836},
                    {"exon_number": 8, "start": 7673534, "end": 7673608},
                    {"exon_number": 9, "start": 7670608, "end": 7670714},
                    {"exon_number": 10, "start": 7669608, "end": 7669689},
                    {"exon_number": 11, "start": 7668421, "end": 7669689},
                ]
            }
        ]
    },
    "BRCA1": {
        "gene_id": "ENSG00000012048",
        "symbol": "BRCA1",
        "chromosome": "chr17",
        "start": 43044295,
        "end": 43170245,
        "strand": "-",
        "biotype": "protein_coding",
        "description": "BRCA1 DNA repair associated",
        "aliases": ["BRCC1", "FANCS", "RNF53"],
        "transcripts": []
    },
    "EGFR": {
        "gene_id": "ENSG00000146648",
        "symbol": "EGFR",
        "chromosome": "chr7",
        "start": 55019017,
        "end": 55211628,
        "strand": "+",
        "biotype": "protein_coding",
        "description": "epidermal growth factor receptor",
        "aliases": ["ERBB1", "HER1"],
        "transcripts": []
    },
}


def get_demo_gene(symbol: str) -> Optional[GeneInfo]:
    """Get demo gene data for testing without database."""
    data = DEMO_GENES.get(symbol.upper())
    if not data:
        return None
    
    transcripts = []
    for t_data in data.get("transcripts", []):
        exons = [
            Exon(
                exon_number=e["exon_number"],
                start=e["start"],
                end=e["end"],
                is_constitutive=True,
            )
            for e in t_data.get("exons", [])
        ]
        
        transcripts.append(Transcript(
            transcript_id=t_data["transcript_id"],
            gene_id=data["gene_id"],
            chromosome=data["chromosome"],
            start=data["start"],
            end=data["end"],
            strand=Strand(data["strand"]),
            cds_start=t_data.get("cds_start"),
            cds_end=t_data.get("cds_end"),
            is_mane_select=t_data.get("is_mane_select", False),
            exons=exons,
        ))
    
    return GeneInfo(
        gene_id=data["gene_id"],
        symbol=data["symbol"],
        chromosome=data["chromosome"],
        start=data["start"],
        end=data["end"],
        strand=Strand(data["strand"]),
        biotype=data["biotype"],
        description=data["description"],
        aliases=data.get("aliases", []),
        transcripts=transcripts,
    )

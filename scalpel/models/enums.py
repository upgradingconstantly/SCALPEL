"""
Core enumerations for SCALPEL.
"""

from enum import Enum


class TargetType(str, Enum):
    """How the target is specified."""
    GENE_SYMBOL = "gene_symbol"
    ENSEMBL_ID = "ensembl_id"
    GENOMIC_COORDINATES = "genomic_coordinates"
    SEQUENCE = "sequence"
    TRANSCRIPT_ID = "transcript_id"


class Genome(str, Enum):
    """Supported reference genomes."""
    # Human
    HUMAN_GRCH38 = "GRCh38"
    HUMAN_GRCH37 = "GRCh37"
    # Mouse
    MOUSE_GRCM39 = "GRCm39"
    MOUSE_GRCM38 = "GRCm38"
    # Zebrafish
    ZEBRAFISH_GRCZ11 = "GRCz11"
    # Rat
    RAT_MRATBN7 = "mRatBN7.2"
    # Drosophila (fruit fly)
    DROSOPHILA_BDGP6 = "BDGP6"
    # C. elegans (worm)
    CELEGANS_WBCEL235 = "WBcel235"
    # Pig
    PIG_SSCROFA11 = "Sscrofa11.1"
    # Arabidopsis (plant)
    ARABIDOPSIS_TAIR10 = "TAIR10"

    @classmethod
    def from_string(cls, value: str) -> "Genome":
        """Parse genome from string, handling common aliases."""
        aliases = {
            "hg38": cls.HUMAN_GRCH38,
            "hg19": cls.HUMAN_GRCH37,
            "mm39": cls.MOUSE_GRCM39,
            "mm10": cls.MOUSE_GRCM38,
            "danrer11": cls.ZEBRAFISH_GRCZ11,
            "dm6": cls.DROSOPHILA_BDGP6,
            "ce11": cls.CELEGANS_WBCEL235,
            "susscr11": cls.PIG_SSCROFA11,
            "tair10": cls.ARABIDOPSIS_TAIR10,
        }
        # Try direct match first
        for genome in cls:
            if genome.value.lower() == value.lower():
                return genome
        # Try aliases
        if value.lower() in aliases:
            return aliases[value.lower()]
        raise ValueError(f"Unknown genome: {value}")


class EditModality(str, Enum):
    """CRISPR editing modality."""
    KNOCKOUT = "knockout"              # SpCas9 DSB-mediated
    INTERFERENCE = "interference"      # dCas9-KRAB (CRISPRi)
    ACTIVATION = "activation"          # dCas9-VP64/VPR (CRISPRa)
    BASE_EDIT_CBE = "base_edit_cbe"    # Cytosine base editor (C→T)
    BASE_EDIT_ABE = "base_edit_abe"    # Adenine base editor (A→G)
    BASE_EDIT_GBE = "base_edit_gbe"    # Glycosylase base editor (C→G)
    PRIME_EDIT = "prime_edit"          # Prime editing
    CRISPROFF = "crisproff"            # Epigenetic silencing (DNA methylation)
    DUAL_NICKASE = "dual_nickase"      # Paired nickases for high specificity


class CasVariantType(str, Enum):
    """Supported Cas protein variants."""
    SPCAS9 = "SpCas9"           # S. pyogenes Cas9 (NGG)
    SPCAS9_NG = "SpCas9-NG"     # SpCas9 variant (NG)
    SPCAS9_VQR = "SpCas9-VQR"   # SpCas9 variant (NGA)
    SPRY = "SpRY"               # SpCas9 variant (NRN/NYN - almost PAM-less)
    XCAS9 = "xCas9"             # Expanded PAM (NG, GAA, GAT)
    SACAS9 = "SaCas9"           # S. aureus Cas9 (NNGRRT)
    SACAS9_KKH = "SaCas9-KKH"   # Relaxed SaCas9 (NNNRRT)
    CJCAS9 = "CjCas9"           # C. jejuni Cas9 (NNNNRYAC) - smallest
    CAS12A = "Cas12a"           # Cpf1 (TTTV)
    CAS12A_RR = "Cas12a-RR"     # Cas12a variant (TYCV)


class PAMPosition(str, Enum):
    """Position of PAM relative to spacer."""
    THREE_PRIME = "3prime"  # PAM after spacer (Cas9)
    FIVE_PRIME = "5prime"   # PAM before spacer (Cas12a)


class Strand(str, Enum):
    """DNA strand."""
    PLUS = "+"
    MINUS = "-"


class RedFlagSeverity(str, Enum):
    """Severity levels for red flags."""
    CRITICAL = "critical"
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class ExperimentType(str, Enum):
    """Types of CRISPR experiments."""
    IN_VITRO = "in_vitro"
    CELL_LINE = "cell_line"
    PRIMARY_CELLS = "primary_cells"
    ORGANOID = "organoid"
    IN_VIVO = "in_vivo"


class ValidationTier(str, Enum):
    """Validation assay tiers."""
    TIER1_GENOTYPING = "genotyping"
    TIER2_MOLECULAR = "molecular"
    TIER3_FUNCTIONAL = "functional"

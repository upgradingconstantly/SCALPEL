# SCALPEL

**S**mart **C**RISPR **A**nalysis **L**ibrary for **P**recision **E**diting and **L**ocus-targeting

A comprehensive computational CRISPR design platform for gRNA design, off-target analysis, experiment planning, and red flag detection across multiple editing modalities.

## Features

### Multi-Modality Support
- **CRISPR Knockout**: SpCas9-mediated double-strand break for gene disruption
- **CRISPRi (Interference)**: dCas9-KRAB transcriptional repression
- **CRISPRa (Activation)**: dCas9-VP64/VPR transcriptional activation
- **Base Editing (CBE)**: Cytosine base editor (C→T conversion)
- **Base Editing (ABE)**: Adenine base editor (A→G conversion)
- **Prime Editing**: Precise insertions, deletions, and substitutions

### Cas Protein Variants
| Variant | PAM | Description |
|---------|-----|-------------|
| SpCas9 | NGG | Standard *S. pyogenes* Cas9 |
| SpCas9-NG | NG | Relaxed PAM variant |
| SpCas9-VQR | NGA | Alternative PAM variant |
| SaCas9 | NNGRRT | *S. aureus* Cas9 (smaller) |
| Cas12a | TTTV | Cpf1, staggered cuts |
| Cas12a-RR | TYCV | Cas12a PAM variant |

### Efficiency Scoring
- **Ensemble scoring**: Combines rule-based (Doench) and sequence feature analysis
- **Position-specific scoring**: Evaluates nucleotide preferences at each position
- **Modality-specific adjustments**: Optimized scoring for each editing approach
- **Confidence intervals**: Uncertainty quantification for predictions

### Off-Target Analysis
- **CFD scoring**: Cutting Frequency Determination algorithm
- **Configurable mismatch tolerance**: Search with 0-4 mismatches
- **Risk categorization**: Sites ranked by cutting probability
- **Specificity scoring**: Overall guide specificity metric

### Red Flag Detection
Automatic detection of potential issues:
- **Sequence issues**: Homopolymers, extreme GC content, Pol III terminators
- **Genomic context**: Repeat elements, segmental duplications
- **Off-target concerns**: High-risk sites, exonic off-targets
- **Gene-level issues**: Essential genes, paralogs
- **SNP overlaps**: Variants that may affect targeting

### Experiment Planning
- **Positive controls**: Validated control guides for each modality
- **Negative controls**: Non-targeting scramble sequences
- **Tiered validation assays**: Genotyping → Molecular → Functional
- **Troubleshooting guides**: Common failure modes and solutions
- **Provenance tracking**: Reproducibility via input hashing

## Supported Genomes

| Species | Genome | Aliases |
|---------|--------|---------|
| Human | GRCh38 | hg38 |
| Human | GRCh37 | hg19 |
| Mouse | GRCm39 | mm39 |
| Mouse | GRCm38 | mm10 |
| Zebrafish | GRCz11 | danrer11 |
| Rat | mRatBN7.2 | - |

## Installation

### Requirements
- Python 3.10+
- Dependencies managed via pip

### Install from source
```bash
git clone https://github.com/upgradingconstantly/SCALPEL.git
cd SCALPEL
pip install -e .
```

### Install with development dependencies
```bash
pip install -e ".[dev]"
```

## Quick Start

### CLI Usage

#### Design gRNAs for a gene
```bash
scalpel design --gene TP53 --modality knockout --genome GRCh38
```

#### Design with specific Cas variant
```bash
scalpel design --gene BRCA1 --cas SpCas9-NG --n-guides 20
```

#### Target by genomic coordinates
```bash
scalpel design --coordinates chr17:7668421-7687490 --genome GRCh38
```

#### Analyze off-targets
```bash
scalpel offtarget --spacer ATCGATCGATCGATCGATCG --genome GRCh38 --max-mismatches 3
```

#### Generate experiment plan
```bash
scalpel plan --input design_results.json --output experiment_plan.md
```

#### Pipe-friendly workflows
```bash
# Full pipeline
scalpel design --gene BRCA1 | scalpel offtarget | scalpel plan

# Save intermediate results
scalpel design --gene MYC --output myc_guides.json
scalpel offtarget --input myc_guides.json --output myc_offtargets.json
```

#### Output formats
```bash
# JSON (default, pipe-friendly)
scalpel design --gene TP53 --format json

# Tab-separated values
scalpel design --gene TP53 --format tsv

# Rich table (terminal display)
scalpel design --gene TP53 --format table

# Streaming JSONL (one guide per line)
scalpel design --gene TP53 --stream
```

#### View configuration and available options
```bash
scalpel info
```

### Python API

```python
from scalpel.design import SpacerExtractor
from scalpel.design.efficiency import EnsembleScorer
from scalpel.offtarget import OffTargetSearcher, RiskCalculator
from scalpel.planning import ExperimentPlanner
from scalpel.core.red_flags import detect_red_flags, summarize_red_flags
from scalpel.models.enums import Genome, EditModality, CasVariantType

# Extract spacers from a sequence
extractor = SpacerExtractor(CasVariantType.SPCAS9)
spacers = extractor.extract_spacers(
    sequence="ATCG...",  # Your target sequence
    chromosome="chr17",
    start_position=7668421
)

# Score spacers for efficiency
scorer = EnsembleScorer()
scored_guides = scorer.score_batch(
    spacers,
    modality="knockout",
    n_top=10  # Return top 10 guides
)

# Analyze off-targets
searcher = OffTargetSearcher(Genome.HUMAN_GRCH38, max_mismatches=4)
for guide in scored_guides:
    analysis = searcher.search(guide.spacer.spacer_sequence)
    risk_calc = RiskCalculator(Genome.HUMAN_GRCH38)
    summary = risk_calc.generate_summary(analysis)
    print(f"{guide.spacer.spacer_sequence}: {summary['specificity_score']}% specificity")

# Detect red flags
for guide in scored_guides:
    flags = detect_red_flags(guide.spacer)
    flag_summary = summarize_red_flags(flags)
    if flag_summary['total_flags'] > 0:
        print(f"Warning: {flag_summary['interpretation']}")

# Generate experiment plan
from scalpel.models.data_classes import DesignedGuide
planner = ExperimentPlanner(EditModality.KNOCKOUT)
plan = planner.generate_plan(scored_guides, gene_info=None)
markdown_plan = planner.to_markdown(plan)
print(markdown_plan)
```

## Web Interface

SCALPEL includes a built-in web server with a modern UI.

### Start the web server
```bash
python -m scalpel.api.server
```

Then open http://localhost:8000 in your browser.

### REST API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/health` | GET | Health check |
| `/api/design` | POST | Design gRNAs |
| `/api/design/stream` | POST | Design with SSE streaming |
| `/api/offtarget` | POST | Off-target analysis |
| `/api/plan` | POST | Generate experiment plan |

### Example API request
```bash
curl -X POST http://localhost:8000/api/design \
  -H "Content-Type: application/json" \
  -d '{"gene": "TP53", "modality": "knockout", "genome": "GRCh38", "n_guides": 10}'
```

## Configuration

SCALPEL uses environment variables and a config file for settings:

```bash
# Data directories
export SCALPEL_DATA_DIR=~/.scalpel
export SCALPEL_GENOMES_DIR=~/.scalpel/genomes
export SCALPEL_MODELS_DIR=~/.scalpel/models
```

## Output Format

### Design Output (JSON)
```json
{
  "status": "success",
  "target": {
    "gene": "TP53",
    "chromosome": "chr17",
    "start": 7668421,
    "end": 7687490,
    "genome": "GRCh38",
    "modality": "knockout"
  },
  "guides": [
    {
      "rank": 1,
      "spacer_sequence": "ATCGATCGATCGATCGATCG",
      "pam_sequence": "AGG",
      "strand": "+",
      "genomic_start": 7668500,
      "genomic_end": 7668520,
      "cut_site": 7668517,
      "efficiency_score": 0.85,
      "efficiency_interpretation": "Excellent efficiency",
      "red_flags": {
        "total_flags": 0,
        "interpretation": "No issues detected"
      }
    }
  ]
}
```

## Development

### Run tests
```bash
pytest tests/
```

### Run specific tests
```bash
pytest tests/test_efficiency.py -v
pytest tests/test_modalities.py -v
```

### Code formatting
```bash
ruff check .
ruff format .
```

### Type checking
```bash
mypy scalpel/
```

## Architecture

```
scalpel/
├── api/              # FastAPI web server
├── benchmark/        # Validation and benchmarking
├── core/             # Core utilities (cache, plugins, red flags)
├── design/           # Spacer extraction and efficiency scoring
│   └── efficiency/   # Ensemble scoring system
├── genome/           # Reference genome handling
├── modalities/       # Modality-specific designers
├── models/           # Data classes and enums
├── offtarget/        # Off-target analysis
├── planning/         # Experiment planning
└── cli.py            # Command-line interface
```

## Citation

If you use SCALPEL in your research, please cite:

```
SCALPEL: Smart CRISPR Analysis Library for Precision Editing and Locus-targeting
https://github.com/upgradingconstantly/SCALPEL
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please read our contributing guidelines and submit pull requests to the main repository.

## Support

- Issues: https://github.com/upgradingconstantly/SCALPEL/issues
- Discussions: https://github.com/upgradingconstantly/SCALPEL/discussions

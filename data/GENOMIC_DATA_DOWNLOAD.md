# SCALPEL Genomic Data Requirements

This folder contains instructions for downloading all genomic data required to run SCALPEL in production.

> [!IMPORTANT]
> SCALPEL requires real genomic reference data to function. The files are large (2-3GB each) and must be downloaded separately.

---

## Quick Start

1. Create the genomes directory: `~/.scalpel/genomes/`
2. Download the genome(s) you need from the links below
3. Decompress and place in the appropriate subdirectory

---

## Required Data Files

### Human Genome - GRCh38 (hg38)

**Recommended for most human CRISPR work.**

| File | Size | Source |
|------|------|--------|
| Genome FASTA | ~3 GB | [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) |
| Gene Annotations | ~50 MB | [GENCODE v44](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz) |

**Alternative Ensembl source:**
- [Ensembl GRCh38 FASTA](https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)

**Setup:**
```bash
mkdir -p ~/.scalpel/genomes/grch38
cd ~/.scalpel/genomes/grch38

# Download genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Download annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

# Create index (requires samtools)
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

---

### Human Genome - GRCh37 (hg19)

**Use for legacy data or compatibility with older databases.**

| File | Size | Source |
|------|------|--------|
| Genome FASTA | ~3 GB | [Ensembl GRCh37](https://ftp.ensembl.org/pub/grch37/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz) |
| Gene Annotations | ~40 MB | [GENCODE v44lift37](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz) |

**Setup:**
```bash
mkdir -p ~/.scalpel/genomes/grch37
cd ~/.scalpel/genomes/grch37

# Download genome
wget https://ftp.ensembl.org/pub/grch37/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# Download annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz
gunzip gencode.v44lift37.annotation.gtf.gz

# Create index
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
```

---

### Mouse Genome - GRCm39 (mm39)

**For mouse CRISPR experiments.**

| File | Size | Source |
|------|------|--------|
| Genome FASTA | ~2.7 GB | [Ensembl GRCm39](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz) |
| Gene Annotations | ~35 MB | [GENCODE vM33](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz) |

**Setup:**
```bash
mkdir -p ~/.scalpel/genomes/grcm39
cd ~/.scalpel/genomes/grcm39

# Download genome
wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Download annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz
gunzip gencode.vM33.annotation.gtf.gz

# Create index
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
```

---

## Building the Gene Database

After downloading annotations, build the gene lookup database:

```bash
# From SCALPEL project root
python -m scalpel.genome.build_gene_db --genome grch38 --gtf ~/.scalpel/genomes/grch38/gencode.v44.annotation.gtf
python -m scalpel.genome.build_gene_db --genome grch37 --gtf ~/.scalpel/genomes/grch37/gencode.v44lift37.annotation.gtf
python -m scalpel.genome.build_gene_db --genome grcm39 --gtf ~/.scalpel/genomes/grcm39/gencode.vM33.annotation.gtf
```

---

## Optional: Off-Target Index

For faster off-target searching, pre-compute the PAM site index:

```bash
# This may take 1-2 hours per genome
python -m scalpel.offtarget.build_index --genome grch38
python -m scalpel.offtarget.build_index --genome grch37
python -m scalpel.offtarget.build_index --genome grcm39
```

---

## Directory Structure

After setup, your genomes directory should look like:

```
~/.scalpel/genomes/
├── grch38/
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
│   ├── gencode.v44.annotation.gtf
│   ├── genes.db
│   └── offtargets.duckdb (optional)
├── grch37/
│   ├── Homo_sapiens.GRCh37.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai
│   ├── gencode.v44lift37.annotation.gtf
│   ├── genes.db
│   └── offtargets.duckdb (optional)
└── grcm39/
    ├── Mus_musculus.GRCm39.dna.primary_assembly.fa
    ├── Mus_musculus.GRCm39.dna.primary_assembly.fa.fai
    ├── gencode.vM33.annotation.gtf
    ├── genes.db
    └── offtargets.duckdb (optional)
```

---

## Verification

Test your setup:

```bash
scalpel design --gene TP53 --genome grch38 --n-guides 5
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Genome FASTA not found" | Check file is in correct directory with correct naming |
| "Gene not found" | Ensure gene database is built from annotations |
| Slow off-target search | Build the pre-computed off-target index |
| Permission errors | Ensure read permissions on genome files |

---

## Additional Resources

- [Ensembl FTP](https://ftp.ensembl.org/pub/)
- [GENCODE](https://www.gencodegenes.org/)
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
- [UCSC Genome Browser Downloads](https://hgdownload.soe.ucsc.edu/downloads.html)

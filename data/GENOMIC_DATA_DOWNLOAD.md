# SCALPEL Genomic Data Requirements

This folder contains instructions for downloading all genomic data required to run SCALPEL in production.

> [!IMPORTANT]
> SCALPEL requires real genomic reference data to function. The files are large (2-3GB each) and must be downloaded separately from their official sources.

---

## Supported Genomes

SCALPEL supports the following reference genomes:

| Genome | Species | Alias | Status |
|--------|---------|-------|--------|
| **GRCh38** | Human | hg38 | ✅ Current |
| **GRCh37** | Human | hg19 | ✅ Legacy |
| **GRCm39** | Mouse | mm39 | ✅ Current |
| **GRCm38** | Mouse | mm10 | ✅ Legacy |
| **GRCz11** | Zebrafish | danRer11 | ✅ Current |
| **mRatBN7.2** | Rat | rn7 | ✅ Current |

---

## Quick Start

```bash
# 1. Create the genomes directory
mkdir -p ~/.scalpel/genomes

# 2. Download your required genome(s) using the instructions below

# 3. Build gene database (after downloading)
python -m scalpel.genome.build_gene_db --genome grch38

# 4. Test your setup
scalpel design --gene TP53 --genome grch38 --n-guides 5
```

---

## Human Genomes

### GRCh38 (hg38) - Recommended

**Source:** NCBI + GENCODE v44

| File | Download Link |
|------|---------------|
| Genome FASTA (~833 MB compressed) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) |
| Gene Annotations GTF | [Download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/grch38
cd ~/.scalpel/genomes/grch38

# Genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

# Index (requires samtools)
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```
</details>

---

### GRCh37 (hg19) - Legacy

**Source:** Ensembl Release 112 + GENCODE v44lift37

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/grch37/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/grch37
cd ~/.scalpel/genomes/grch37

# Genome
wget https://ftp.ensembl.org/pub/grch37/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# Annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gtf.gz
gunzip gencode.v44lift37.annotation.gtf.gz

# Index
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
```
</details>

---

## Mouse Genomes

### GRCm39 (mm39) - Recommended

**Source:** Ensembl Release 112 + GENCODE M33

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/grcm39
cd ~/.scalpel/genomes/grcm39

# Genome
wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz
gunzip gencode.vM33.annotation.gtf.gz

# Index
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
```
</details>

---

### GRCm38 (mm10) - Legacy

**Source:** Ensembl Release 102 + GENCODE M25

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/grcm38
cd ~/.scalpel/genomes/grcm38

# Genome
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz

# Index
samtools faidx Mus_musculus.GRCm38.dna.primary_assembly.fa
```
</details>

---

## Zebrafish Genome

### GRCz11 (danRer11)

**Source:** Ensembl Release 112

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ensembl.org/pub/release-112/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/grcz11
cd ~/.scalpel/genomes/grcz11

# Genome
wget https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz

# Annotations
wget https://ftp.ensembl.org/pub/release-112/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz
gunzip Danio_rerio.GRCz11.112.gtf.gz

# Index
samtools faidx Danio_rerio.GRCz11.dna.primary_assembly.fa
```
</details>

---

## Rat Genome

### mRatBN7.2 (rn7)

**Source:** Ensembl Release 112

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ensembl.org/pub/release-112/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/mratbn7.2
cd ~/.scalpel/genomes/mratbn7.2

# Genome
wget https://ftp.ensembl.org/pub/release-112/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz
gunzip Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz

# Annotations
wget https://ftp.ensembl.org/pub/release-112/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.112.gtf.gz
gunzip Rattus_norvegicus.mRatBN7.2.112.gtf.gz

# Index
samtools faidx Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa
```
</details>

---

## Building Gene Databases

After downloading annotations, build the gene lookup database:

```bash
# Human
python -m scalpel.genome.build_gene_db --genome grch38 --gtf ~/.scalpel/genomes/grch38/gencode.v44.annotation.gtf
python -m scalpel.genome.build_gene_db --genome grch37 --gtf ~/.scalpel/genomes/grch37/gencode.v44lift37.annotation.gtf

# Mouse  
python -m scalpel.genome.build_gene_db --genome grcm39 --gtf ~/.scalpel/genomes/grcm39/gencode.vM33.annotation.gtf
python -m scalpel.genome.build_gene_db --genome grcm38 --gtf ~/.scalpel/genomes/grcm38/gencode.vM25.annotation.gtf

# Zebrafish
python -m scalpel.genome.build_gene_db --genome grcz11 --gtf ~/.scalpel/genomes/grcz11/Danio_rerio.GRCz11.112.gtf

# Rat
python -m scalpel.genome.build_gene_db --genome mratbn7.2 --gtf ~/.scalpel/genomes/mratbn7.2/Rattus_norvegicus.mRatBN7.2.112.gtf
```

---

## Optional: Off-Target Index

For faster off-target searching, pre-compute the PAM site index:

```bash
# This takes 1-2 hours per genome
python -m scalpel.offtarget.build_index --genome grch38
python -m scalpel.offtarget.build_index --genome grcm39
```

---

## Expected Directory Structure

```
~/.scalpel/genomes/
├── grch38/
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
│   ├── gencode.v44.annotation.gtf
│   └── genes.db
├── grch37/
│   └── ...
├── grcm39/
│   └── ...
├── grcm38/
│   └── ...
├── grcz11/
│   └── ...
└── mratbn7.2/
    └── ...
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Genome FASTA not found" | Verify file is in `~/.scalpel/genomes/{genome}/` with correct naming |
| "Gene not found" | Build gene database from annotations using `build_gene_db` |
| Slow off-target search | Build pre-computed index using `build_index` |
| Download fails | Check network; try alternative mirror (Ensembl ↔ NCBI) |

---

## Additional Resources

| Resource | URL |
|----------|-----|
| Ensembl FTP | https://ftp.ensembl.org/pub/ |
| GENCODE | https://www.gencodegenes.org/human/ |
| NCBI GenBank | https://ftp.ncbi.nlm.nih.gov/genomes/ |
| UCSC Genome Browser | https://hgdownload.soe.ucsc.edu/downloads.html |

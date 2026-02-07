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
| **BDGP6** | Drosophila | dm6 | ✅ Current |
| **WBcel235** | C. elegans | ce11 | ✅ Current |
| **Sscrofa11.1** | Pig | susScr11 | ✅ Current |
| **TAIR10** | Arabidopsis | - | ✅ Current |

---

## Quick Start

```bash
# 1. Create the genomes directory
mkdir -p ~/.scalpel/genomes

# 2. Download your required genome(s) using the instructions below

# 3. (Optional but recommended for off-target speed) build off-target index
python -m scalpel.offtarget.build_index --genome GRCh38

# 4. Test your setup
scalpel design --gene TP53 --genome GRCh38 --n-guides 5
```

Notes:
- Baseline guide scoring and off-target lookup do not require model training.
- `build_gene_db` is not currently shipped as a supported command in this codebase.

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

## Drosophila Genome

### BDGP6 (dm6)

**Source:** Ensembl Release 112

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-112/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ensembl.org/pub/release-112/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.112.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/bdgp6
cd ~/.scalpel/genomes/bdgp6

# Genome
wget https://ftp.ensembl.org/pub/release-112/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz

# Annotations
wget https://ftp.ensembl.org/pub/release-112/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.112.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.46.112.gtf.gz

# Index
samtools faidx Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
```
</details>

---

## C. elegans Genome

### WBcel235 (ce11)

**Source:** Ensembl Release 112

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-112/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ensembl.org/pub/release-112/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.112.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/wbcel235
cd ~/.scalpel/genomes/wbcel235

# Genome
wget https://ftp.ensembl.org/pub/release-112/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz

# Annotations
wget https://ftp.ensembl.org/pub/release-112/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.112.gtf.gz
gunzip Caenorhabditis_elegans.WBcel235.112.gtf.gz

# Index
samtools faidx Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
```
</details>

---

## Pig Genome

### Sscrofa11.1

**Source:** Ensembl Release 112

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensembl.org/pub/release-112/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ensembl.org/pub/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/sscrofa11.1
cd ~/.scalpel/genomes/sscrofa11.1

# Genome
wget https://ftp.ensembl.org/pub/release-112/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz

# Annotations
wget https://ftp.ensembl.org/pub/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz
gunzip Sus_scrofa.Sscrofa11.1.112.gtf.gz

# Index
samtools faidx Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
```
</details>

---

## Arabidopsis Genome

### TAIR10

**Source:** Ensembl Plants Release 59

| File | Download Link |
|------|---------------|
| Genome FASTA | [Download](https://ftp.ensemblgenomes.org/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz) |
| Gene Annotations GTF | [Download](https://ftp.ensemblgenomes.org/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz) |

<details>
<summary><b>Setup Commands</b></summary>

```bash
mkdir -p ~/.scalpel/genomes/tair10
cd ~/.scalpel/genomes/tair10

# Genome
wget https://ftp.ensemblgenomes.org/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Annotations
wget https://ftp.ensemblgenomes.org/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.59.gtf.gz

# Index
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```
</details>

---

## Gene Database Status

`build_gene_db` command references are currently unsupported in this repository revision.
Keep GTF files for future use, but do not rely on `python -m scalpel.genome.build_gene_db ...` at this time.

---

## Optional: Off-Target Index

For faster off-target searching, pre-compute the PAM site index:

```bash
# This takes 1-2 hours per genome
python -m scalpel.offtarget.build_index --genome GRCh38
python -m scalpel.offtarget.build_index --genome GRCm39
```

Runtime lookup order for off-target DB:
1. `SCALPEL_OFFTARGET_DB` (explicit path)
2. `~/.scalpel/offtargets/{genome}.duckdb`
3. `~/.scalpel/genomes/{genome}/offtargets.duckdb`

---

## Expected Directory Structure

```
~/.scalpel/offtargets/
├── GRCh38.duckdb
└── GRCm39.duckdb

~/.scalpel/genomes/
├── grch38/
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
│   ├── gencode.v44.annotation.gtf
│   └── offtargets.duckdb   # optional fallback location
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
| "Gene not found" | Verify symbol/genome pairing and ensure target annotations are available for that genome |
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

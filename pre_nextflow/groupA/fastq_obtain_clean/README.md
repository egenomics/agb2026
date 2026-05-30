# fastq_obtain_clean

This folder contains the scripts used to **retrieve SRA accession numbers, dowload raw sequencing data, and convert it to FASTQ format** from the NCBI SRA database.

The three scripts run sequenctially and take as input the celaned metadata file produced by `cleaning_metadata.py`.

## Files 

| File | Description |
|------|-------------|
| `retrieve_sra.sh` | Queries NCBI SRA to resolve a run accession (SRR/ERR/DRR) for each sample and adds it to the metadata as a `sample-id` column. |
| `download_sra.sh` | Downloads raw `.sra` files from NCBI using `prefetch` for each accession in the `sample-id` column. |
| `fasterq_dump.sh` | Converts `.sra` files to FASTQ format using `fasterq-dump`, handling both single-end and paired-end layouts automatically. |

## Requirements

```bash
# For retrieve_sra.sh
conda create -n retrieve_sra_env -y
conda install -n retrieve_sra_env bioconda::entrez-direct -y

# For download_sra.sh and fasterq_dump.sh 
conda create -n fasterq_dump_env -y
conda install -n fasterq_dump_env bioconda::sra-tools -y
```
> `bioconda::entrez-direct` is required to query NCBI database for search by `sample_name_id column` and then search each SRA code with `esearch` and retrieve the accession number with `efetch`.
> `bioconda::sra-tools` is required to work with `.sra` files to download them with `prefetch` and convert them into `.fastq`files with `fasterq-dump`.

## Pipeline

#### 1) Retriev SRA accessions (`retrieve_sra.sh`)

Reads the `sample_name_id` column from `sample_information_cleaned_1834.tsv` and queries NCBI SRA via `esearch` + `efetch` to find the corresponding run accession (SRR/ERR/DRR). The accession is added as a new `sample-id` column. 

Samples that cannot be resolved are flagged as `NA`, reported to the terminal, logged to a separated file and removed from the output.

**Input:** The resulting file from `cleaning_metadata.py`
- `<input>.tsv` # sample_information_cleaned_1834.tsv

**Output:** 
- `<input>_clean.tsv` # sample_information_cleaned_1834_clean.tsv
- `<output>_failed_samples.txt` 

| Flag | Required | Description |
|------|----------|-------------|
| `-i` | Yes | Path to the input metadata TSV |
| `-o` | No | Path to the output file (defaults to `<input>_clean.tsv`) |


#### 2) Download `.sra` files (`dowlnoad_sra.sh`)

Reads the `sample-id` column from `<input>_clean.tsv` (sample_information_cleaned_1834_clean.tsv). Files are dowloaded to a temporary subdirectory, then moved to a flat output directory. Failed dowloads are reported per sampl without interrupting the overall run.

**Input:** The resulting file from `retrieve_sra.sh`
-  `<input>_clean.tsv` # sample_information_cleaned_1834_clean.tsv
- Output directory path

**Output:** Generate a directory with all `.sra` files
- `<OUTPUT_DIR>/<sample_id>.sra` 

| Flag | Required | Description |
|------|----------|-------------|
| `-s` | Yes | Path to the metadata TSV |
| `-o` | Yes | Path to the output directory for `.sra` files |

#### 3) Convert `.sra` files into `.fastq` files (`fasterq_dump.sh`)

Reads sample accessions from `sample-id` column of the metadata, locates the corresponding `.sra` file in the input directory and runs `fasterq-dump --split-3`. The `--split-3` flag automatically handles both single-end (SE) and paired-end (PE) layouts: SE data produces one `.fastq` file, PE data produces two files (`_1.fastq`and `_2.fastq`). Samples whose `.sra` file is not found are skipped with a warning.

**Input:**
- Directory containing `.sra` files (output of `download_sra.sh`).
- Base output directory for FASTQ files.
- Metadata TSV with a `sample-id` column: `sample_information_cleaned_1834_clean.tsv`
**Output:**
- `<OUTPUT_DIR>/inter_data/seqs/splitted/` — FASTQ files per sample.

| Flag | Required | Description |
|------|----------|-------------|
| `-i` | Yes | Directory containing `.sra` files |
| `-o` | Yes | Base output directory for FASTQ files |
| `-s` | Yes | Path to the metadata TSV |

### Usage

``bash
# 1) Retrieve SRA accessions
conda activate retireve_sra_env
bash retrieve_sra.sh -i sample_information_cleaned_1834.tsv [-o output_file.tsv]
# 2) Download .sra files
conda activate fasterq_dump_env
bash download_sra.sh -s sample_information_cleaned_1834_clean.tsv -o data/sra
# 3) Convert .sra into .fastq files
conda activate fasterq_dump_env
bash fasterq_dump.sh -i data/sra -o data/ -s sample_information_cleaned_1834_clean.tsv
```
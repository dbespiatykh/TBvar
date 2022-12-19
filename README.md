<img align ="left" src=logo.svg width=250px style="padding-right: 25px; padding-top: 25px;">

# A snakemake workflow for variant calling and lineage barcoding of the _Mycobacterium tuberculosis_ samples.

- [Installation](#installation)
- [Usage](#usage)

## Installation

Use the [Conda](https://docs.conda.io/en/latest/) package manager and [BioConda](https://bioconda.github.io/index.html) channel to install TBvar.

If you do not have conda installed do the following:

```bash
# Download Conda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Set permissions
chmod -X Miniconda3-latest-Linux-x86_64.sh
# Install
bash Miniconda3-latest-Linux-x86_64.sh
```

Set up channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Get **TBvar** pipeline:

```bash
git clone https://github.com/dbespiatykh/TBvar.git
```

Install required dependencies:

```bash
conda install -c conda-forge mamba
mamba env create --file environment.yml
```

## Usage

Activate **TBvar** environment:

```bash
conda activate TBvar
cd TBvar
```

:point_right: In `config` folder edit `config.yml` and add your `samples.tsv` table location, it should be formatted like this:

| Run_accession | R1                             | R2                             |
| ------------- | ------------------------------ | ------------------------------ |
| SRR2024996    | /path/to/SRR2024996_1.fastq.gz | /path/to/SRR2024996_2.fastq.gz |
| SRR2024925    | /path/to/SRR2024925_1.fastq.gz | /path/to/SRR2024925_2.fastq.gz |
| SRR12882189   | /path/to/SRR12882189.fastq.gz  |                                |

**Run_accession** - Run accession number or sample name;\
**R1** - Path to the first read pair;\
**R2** - Path to the second read pair.

Run pipeline:

```bash
snakemake --conda-frontend mamba --use-conda -j 48 -c 48 --max-threads 48 -k --rerun-incomplete
```

It is recommended to use dry run if you are running pipeline for the first time, to see if everything is in working order, for this you can use `-n` flag:

```bash
snakemake --conda-frontend mamba -np
```

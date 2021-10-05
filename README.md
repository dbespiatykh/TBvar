- [Description](#Description)
- [Installation](#Installation)
- [Usage](#Usage)

## Description

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

Get [SpolLineages](https://github.com/dcouvin/SpolLineages):

```bash
cd TBvar && mkdir bin && cd $_
git clone https://github.com/dcouvin/SpolLineages.git
cd SpolLineages/Binary_Mask2
gcc Mask2.c -o Mask2
```
Install required dependencies:

```bash
cd TBvar
conda env create --file environment.yml
```

## Usage

Activate **TBvar** environment:

```bash
conda activate TBvar
cd TBvar
```

:point_right: In project folder make `FASTQ` folder:

```bash
mkdir FASTQ
```

and move to it (symlink) your paired- or single-end [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) reads. Paired-end FASTQ suffix should be: `_1.fastq.gz` and `_2.fastq.gz`, for single-end it should be: `_1.fastq.gz`.

> If you are using both PE and SE reads, firts time you run the pipleine it will throw an error `AmbiguousRuleException:...`, just ignore it and run it again.

Project folder should have the following structure:

```bash
📂TBvar/
├── FASTQ
│   ├── SampleA_1.fastq.gz
|   ├── SampleA_2.fastq.gz
|   ├── SampleB_1.fastq.gz
|   └── ...
├── LICENSE
├── README.md
├── Snakefile
├── bin
│   └── SpolLineages
│       ├── Binary_Mask2
|       ├── ...
|       ...
├── environment.yml
├── envs
│   ├── barcoding.yaml
│   ├── bcftools.yaml
│   ├── rtg_tools.yaml
│   ├── snpit.yaml
│   ├── spollineages.yaml
│   ├── spotyping.yaml
│   └── tb_profiler.yaml
├── ref
│   ├── NC_000962.3.dict
│   ├── NC_000962.3.fna
│   ├── NC_000962.3.fna.amb
│   ├── NC_000962.3.fna.ann
│   ├── NC_000962.3.fna.bwt
│   ├── NC_000962.3.fna.fai
│   ├── NC_000962.3.fna.pac
│   ├── NC_000962.3.fna.sa
│   └── tbdb
│       ├── tbdb.ann.txt
│       ├── tbdb.barcode.bed
│       ├── tbdb.bed
│       ├── tbdb.confidence.csv
│       ├── tbdb.csv
│       ├── tbdb.dr.json
│       ├── tbdb.fasta
│       ├── tbdb.fasta.fai
│       ├── tbdb.gff
│       ├── tbdb.version.json
│       └── tbdb.watchlist.csv
└──scripts
    ├── aggregate.py
    ├── concat.py
    ├── gather_spol.py
    ├── read_vcf_v3.py
    └── rtg_concat.py
```

Run pipeline:

```bash
snakemake --conda-frontend mamba --use-conda -j 48 -c 48 --max-threads 48 -k --rerun-incomplete
```

It is recommended to use dry run if you are running pipeline for the first time, to see if everything is in working order, for this you can use `-n` flag:

```bash
snakemake --conda-frontend mamba --use-conda -j 48 -c 48 --max-threads 48 -k --rerun-incomplete -np
```

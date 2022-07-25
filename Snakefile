## Configuration file ----------------------------------------

configfile: "config.yml"

## Set samples -----------------------------------------------

(SAMPLES,) = glob_wildcards("FASTQ/{smp}_1.fastq.gz")

## Target rule -----------------------------------------------

rule all:
    input:
        "Results/all.snps.pass.barcoded.tsv",

## Workflow rules --------------------------------------------

include: "rules/reference.smk",
include: "rules/mapping.smk",
include: "rules/calling.smk",
include: "rules/barcoding.smk"

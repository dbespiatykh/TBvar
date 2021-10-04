import os
import glob
import time
import shutil
import ntpath
import datetime
import numpy as np
import pandas as pd

t_start = time.time()

## Setting order of rules
ruleorder: BWA_MEM_pe > BWA_MEM_se
ruleorder: spotyping_pe > spotyping_se


## ------------------------------------------------------------------------------------ ##
## Globals
## ------------------------------------------------------------------------------------ ##

(SAMPLES,) = glob_wildcards("FASTQ/{smp}_1.fastq.gz")


## ------------------------------------------------------------------------------------ ##
## Workflow rule
## ------------------------------------------------------------------------------------ ##
rule all:
    input:
        "VCF/all.snps.pass.vcf.gz",
        "results/qc/multiqc.html",
        "results/mtbc_barcodes.tsv",
        "stats/rtg_vcfstats/vcf_stats.tsv",
        "results/tb_profiler/tbprofiler.txt",
        "results/mtbc_spoligotypes.tsv",
        "results/mtbc_spoligotypes_annotated.tsv",
        "results/tbvar_results.xlsx",

## ------------------------------------------------------------------------------------ ##
## Rules
## ------------------------------------------------------------------------------------ ##


## Map PE reads
rule BWA_MEM_pe:
    input:
        reads=["FASTQ/{smp}_1.fastq.gz", "FASTQ/{smp}_2.fastq.gz"],
        idx="ref/NC_000962.3.fna",
    output:
        temp("BAM/{smp}.raw.bam")
    log:
        "logs/bwa_mem/{smp}.log"
    params:
        index="ref/NC_000962.3.fna",
        extra=r"-R '@RG\tID:{smp}\tSM:{smp}'",
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.74.0/bio/bwa/mem"


## Map SE reads
rule BWA_MEM_se:
    input:
        reads="FASTQ/{smp}_1.fastq.gz",
        idx="ref/NC_000962.3.fna",
    output:
        temp("BAM/{smp}.raw.bam")
    log:
        "logs/bwa_mem/{smp}.log"
    params:
        index="ref/NC_000962.3.fna",
        extra=r"-R '@RG\tID:{smp}\tSM:{smp}'",
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.74.0/bio/bwa/mem"


## Remove duplicate reads
rule mark_duplicates:
    input:
        "BAM/{smp}.raw.bam"
    output:
        bam="BAM/{smp}.bam",
        metrics="stats/picard/{smp}.metrics.txt"
    log:
        "logs/picard/dedup/{smp}.log"
    params:
        "REMOVE_DUPLICATES=true"
    resources:
        mem_mb=10240
    wrapper:
        "0.77.0/bio/picard/markduplicates"


## Index BAM files
rule samtools_index:
    input:
        "BAM/{smp}.bam",
    output:
        "BAM/{smp}.bam.bai",
    log:
        "logs/samtools/index/{smp}.log",
    wrapper:
        "0.74.0/bio/samtools/index"


## Call variants
rule haplotype_caller:
    input:
        bam="BAM/{smp}.bam",
        idx="BAM/{smp}.bam.bai",
        ref="ref/NC_000962.3.fna",
    output:
        gvcf="VCF/{smp}.g.vcf.gz",
    log:
        "logs/gatk/haplotypecaller/{smp}.log"
    params:
        extra="-ploidy 1 -mbq 20",
    resources:
        mem_mb=10240
    wrapper:
        "0.77.0/bio/gatk/haplotypecaller"


## Import VCFs to GenomicsDB
rule genomics_db_import:
    input:
        gvcfs=expand("VCF/{smp}.g.vcf.gz", smp=SAMPLES),
    output:
        db=directory("db"),
    log:
        "logs/gatk/genomicsdbimport.log"
    params:
        intervals="NC_000962.3",
        db_action="create",
    resources:
        mem_mb=200000
    wrapper:
        "0.77.0/bio/gatk/genomicsdbimport"


## Convert GenomicsDB to VCF
rule genotype_gvcfs:
    input:
        genomicsdb="db",
        ref="ref/NC_000962.3.fna",
    output:
        vcf="VCF/all.vcf.gz",
    log:
        "logs/gatk/genotypegvcfs.log"
    resources:
        mem_mb=200000
    wrapper:
        "0.77.0/bio/gatk/genotypegvcfs"


## Select only SNPs
rule gatk_select:
    input:
        vcf="VCF/all.vcf.gz",
        ref="ref/NC_000962.3.fna",
    output:
        vcf="VCF/all.snps.vcf.gz",
    log:
        "logs/gatk/select/snps.log"
    params:
        extra="--select-type-to-include SNP",
    resources:
        mem_mb=200000
    wrapper:
        "0.77.0/bio/gatk/selectvariants"


## Filter SNPs
rule gatk_filter:
    input:
        vcf="VCF/all.snps.vcf.gz",
        ref="ref/NC_000962.3.fna",
    output:
        vcf="VCF/all.snps.filtered.vcf.gz",
    log:
        "logs/gatk/filter/snps.log"
    params:
        filters={"Failfilter": "QD < 2.0 || DP < 10 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
    resources:
        mem_mb=200000
    wrapper:
        "0.77.0/bio/gatk/variantfiltration"


## Select only PASS variants
rule gatk_select_passed:
    input:
        vcf="VCF/all.snps.filtered.vcf.gz",
        ref="ref/NC_000962.3.fna",
    output:
        vcf="VCF/all.snps.pass.vcf.gz",
    log:
        "logs/gatk/select/snps.pass.log"
    params:
        extra="--exclude-filtered",
    resources:
        mem_mb=200000
    wrapper:
        "0.77.0/bio/gatk/selectvariants"


## Split VCF to samples
rule gatk_select_samples:
    input:
        vcf="VCF/all.snps.pass.vcf.gz",
        ref="ref/NC_000962.3.fna",
    output:
        vcf="VCF/split/{smp}.vcf.gz",
    log:
        "logs/gatk/select/split/{smp}.log"
    params:
        extra="--sample-name {smp}",
    resources:
        mem_mb=4096
    wrapper:
        "0.77.0/bio/gatk/selectvariants"


## Run fastqc
rule fastqc:
    input:
        "BAM/{smp}.bam",
    output:
        html="stats/fastqc/{smp}.html",
        zip="stats/fastqc/{smp}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{smp}.log"
    threads: 8
    wrapper:
        "0.77.0/bio/fastqc"


## Run samtools-stats
rule samtools_stats:
    input:
        bam="BAM/{smp}.bam",
    output:
        "stats/samtools_stats/{smp}.txt"
    log:
        "logs/samtools_stats/{smp}.log"
    wrapper:
        "0.77.0/bio/samtools/stats"


## Run rtg vcfstats
rule rtg_vcfstats:
    input:
        vcf="VCF/all.snps.pass.vcf.gz",
    output:
        "stats/rtg_vcfstats/all.snps.pass.stats.txt",
    conda:
        "envs/rtg_tools.yaml"
    shell:
        "rtg vcfstats {input.vcf} | sed -e 's/:/\\t/' -e 's/ //g' -e '1,4d' \
            | gawk -F '\\t' '{{OFS = '\\t'; print $2}}' > {output}"


## Make vcf stats table
rule vcfstats_table:
    input:
       "stats/rtg_vcfstats/all.snps.pass.stats.txt",
    output:
        "stats/rtg_vcfstats/vcf_stats.tsv",
    conda:
        "envs/barcoding.yaml"
    script:
        "scripts/rtg_concat.py"


## Run multiqc
rule multiqc:
    input:
        expand(
            ["stats/samtools_stats/{smp}.txt",
            "stats/fastqc/{smp}_fastqc.zip",
            "stats/picard/{smp}.metrics.txt"
            ],
            smp=SAMPLES
            )
    output:
        "results/qc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.77.0/bio/multiqc"


## Barcode samples
rule genotyping:
    input:
        "VCF/split/{smp}.vcf.gz",
    output:
        "results/barcoding/{smp}.tsv",
    conda:
        "envs/barcoding.yaml"
    script:
        "scripts/read_vcf_v3.py"


## Concatenate barcodes
rule concatenate_barcodes:
    input:
       expand("results/barcoding/{smp}.tsv", smp=SAMPLES),
    output:
        "results/mtbc_barcodes.tsv",
    conda:
        "envs/barcoding.yaml"
    script:
        "scripts/concat.py"


## Predicting lineages and drug resistance with TBProfiler
rule run_tb_profiler:
    input:
        vcf="VCF/split/{smp}.vcf.gz",
    output:
        "results/tb_profiler/results/{smp}.results.json",
    conda:
        "envs/tb_profiler.yaml"
    log:
        "logs/tb_profiler/{smp}.log",
    shell:
        "tb-profiler vcf_profile \
                                           --external_db ref/tbdb/tbdb \
                                           --dir results/tb_profiler \
                                           {input.vcf} 2> {log}"


## Collating TBProfiler outputs
rule collate_tb_profiler_results:
    input:
        expand("results/tb_profiler/results/{smp}.results.json", smp=SAMPLES),
    output:
        "results/tb_profiler/tbprofiler.txt",
    conda:
        "envs/tb_profiler.yaml"
    log:
        "logs/tb_profiler/collate.log",
    shell:
        "tb-profiler collate \
                                           --dir results/tb_profiler/results \
                                           --full \
                                           --external_db ref/tbdb/tbdb \
                                           --prefix results/tb_profiler/tbprofiler 2> {log}"



## Run SpoTyping on PE
rule spotyping_pe:
    input:
        r1="FASTQ/{smp}_1.fastq.gz",
        r2="FASTQ/{smp}_2.fastq.gz",
    output:
        "results/spotyping/{smp}.SpoTyping",
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.ndb"),
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.not"),
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.ntf"),
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.nto"),
    conda:
        "envs/spotyping.yaml"
    log:
        "logs/spotyping/{smp}.log",
    shell:
        "(SpoTyping.py --noQuery -O results/spotyping -o {wildcards.smp}.SpoTyping {input.r1} {input.r2}) 2> {log}"


## Run SpoTyping on SE
rule spotyping_se:
    input:
        r1="FASTQ/{smp}_1.fastq.gz",
    output:
        "results/spotyping/{smp}.SpoTyping",
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.ndb"),
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.not"),
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.ntf"),
        temp("results/spotyping/{smp}.SpoTyping.SpoTyping.tmp.0.nto"),
    conda:
        "envs/spotyping.yaml"
    log:
        "logs/spotyping/{smp}.log",
    shell:
        "(SpoTyping.py --noQuery -O results/spotyping -o {wildcards.smp}.SpoTyping {input.r1}) 2> {log}"


## Concatenate spoligotypes
rule concatenate_spoligotypes:
    input:
       expand("results/spotyping/{smp}.SpoTyping", smp=SAMPLES),
    output:
        "results/mtbc_spoligotypes.tsv",
    conda:
        "envs/barcoding.yaml"
    script:
        "scripts/gather_spol.py"


## Annotate spoligotypes
rule spollineages:
    input:
        tab="results/mtbc_spoligotypes.tsv",
    output:
        "results/mtbc_spoligotypes_annotated.tsv",
    log:
        "logs/spollineages.log",
    shell:
        "(java -jar bin/SpolLineages/spollineages.jar -delim tab \
                                                                         -i {input.tab} \
                                                                         -o {output} \
                                                                         -D -pDT bin/SpolLineages/Decision_Tree/ \
                                                                         -E -pEA bin/SpolLineages/Binary_Mask2/ \
                                                                         -a) &> {log}"

## Aggregate results
rule aggregate_results:
    input:
       "stats/rtg_vcfstats/vcf_stats.tsv",
       "results/qc/multiqc_data/multiqc_fastqc.txt",
       "results/qc/multiqc_data/multiqc_samtools_stats.txt",
       "results/qc/multiqc_data/multiqc_picard_dups.txt",
       "results/mtbc_barcodes.tsv",
       "results/mtbc_spoligotypes_annotated.tsv",
       "results/tb_profiler/tbprofiler.txt",
    output:
        "results/tbvar_results.xlsx",
    conda:
        "envs/barcoding.yaml"
    script:
        "scripts/aggregate.py"

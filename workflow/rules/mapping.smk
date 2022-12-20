## Map PE/SE reads
rule bwa_mem_mapping:
    input:
        reads=get_fastq,
        idx=multiext("resources/ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("results/BAM/{sample}.srt.bam"),
    log:
        "logs/bwa/mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",
        sort_order="coordinate",
    threads: config["BWA"]["threads"]
    wrapper:
        "v1.21.0/bio/bwa/mem"


## Remove duplicated reads
rule picard_mark_duplicates:
    input:
        bams="results/BAM/{sample}.srt.bam",
    output:
        bam="results/BAM/{sample}.bam",
        metrics="results/stats/picard/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        "REMOVE_DUPLICATES=true",
    resources:
        mem_mb=config["GATK"]["markdup"]["memory"],
    wrapper:
        "v1.21.0/bio/picard/markduplicates"


## Index BAM files
rule samtools_index:
    input:
        "results/BAM/{sample}.bam",
    output:
        "results/BAM/{sample}.bam.bai",
    log:
        "logs/samtools/index/{sample}.log",
    wrapper:
        "v1.21.0/bio/samtools/index"

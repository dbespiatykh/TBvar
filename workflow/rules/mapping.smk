## Map PE reads
rule bwa_mem_mapping:
    input:
        reads=get_fastq,
        idx=multiext("ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("BAM/{sample}.srt.bam"),
    log:
        "logs/bwa/mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",
        sort_order="coordinate",
    threads: config["BWA"]["threads"]
    wrapper:
        "v1.7.1/bio/bwa/mem"


## Remove duplicated reads
rule picard_mark_duplicates:
    input:
        bams="BAM/{sample}.srt.bam",
    output:
        bam="BAM/{sample}.bam",
        metrics="stats/picard/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        "REMOVE_DUPLICATES=true",
    resources:
        mem_mb=config["GATK"]["markdup"]["memory"],
    wrapper:
        "v1.7.1/bio/picard/markduplicates"


## Index BAM files
rule samtools_index:
    input:
        "BAM/{sample}.bam",
    output:
        "BAM/{sample}.bam.bai",
    log:
        "logs/samtools/index/{sample}.log",
    wrapper:
        "v1.7.1/bio/samtools/index"

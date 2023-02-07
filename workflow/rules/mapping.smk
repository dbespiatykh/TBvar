## Map PE/SE reads
rule bwa_mem_mapping:
    input:
        reads=get_fastq,
        idx=multiext(
            "resources/ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
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
        "v1.23.0/bio/bwa/mem"


## Remove duplicated reads
rule sambamba_mark_duplicates:
    input:
        bams="results/BAM/{sample}.srt.bam",
    output:
        bam=protected("results/BAM/{sample}.bam"),
    params:
        extra="-r",
    log:
        "logs/sambamba-markdup/{sample}.log",
    threads: config["SAMBAMBA"]["threads"]
    wrapper:
        "v1.23.0/bio/sambamba/markdup"


## Index BAM files
rule sambamba_index:
    input:
        "results/BAM/{sample}.bam",
    output:
        "results/BAM/{sample}.bam.bai",
    params:
         extra="",
    log:
        "logs/sambamba-index/{sample}.log",
    threads: config["SAMBAMBA"]["threads"]
    wrapper:
        "v1.23.0/bio/sambamba/index"

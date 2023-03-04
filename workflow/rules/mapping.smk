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
        idx=protected("results/BAM/{sample}.bam.bai"),
    params:
        extra="-r",
    log:
        "logs/sambamba-markdup/{sample}.log",
    threads: config["SAMBAMBA"]["threads"]
    wrapper:
        "v1.23.0/bio/sambamba/markdup"

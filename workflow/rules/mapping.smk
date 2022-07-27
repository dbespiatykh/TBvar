## Setting order of rules
ruleorder: BWA_MEM_pe > BWA_MEM_se > mark_duplicates

## Map PE reads
rule BWA_MEM_pe:
    input:
        reads=[get_fastq(sample)["R1"], get_fastq(sample)["R2"]],
        idx=multiext("ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("BAM/{smp}.srt.bam"),
    log:
        "logs/bwa_mem/{smp}.log",
    params:
        extra=r"-R '@RG\tID:{smp}\tSM:{smp}'",
        sorting="samtools",
        sort_order="coordinate",
    threads: config["bwa_threads"]
    wrapper:
        "v1.7.1/bio/bwa/mem"


## Map SE reads
rule BWA_MEM_se:
    input:
        reads=get_fastq(sample)["R1"],
        idx=multiext("ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("BAM/{smp}.srt.bam"),
    log:
        "logs/bwa_mem/{smp}.log",
    params:
        extra=r"-R '@RG\tID:{smp}\tSM:{smp}'",
        sorting="samtools",
        sort_order="coordinate",
    threads: config["bwa_threads"]
    wrapper:
        "v1.7.1/bio/bwa/mem"


## Remove duplicate reads
rule mark_duplicates:
    input:
        bams="BAM/{smp}.srt.bam",
    output:
        bam="BAM/{smp}.bam",
        metrics="stats/picard/{smp}.metrics.txt",
    log:
        "logs/picard/dedup/{smp}.log",
    params:
        "REMOVE_DUPLICATES=true",
    resources:
        mem_mb=config["markdup_mem"],
    wrapper:
        "v1.7.1/bio/picard/markduplicates"


## Index BAM files
rule samtools_index:
    input:
        "BAM/{smp}.bam",
    output:
        "BAM/{smp}.bam.bai",
    log:
        "logs/samtools/index/{smp}.log",
    wrapper:
        "v1.7.1/bio/samtools/index"

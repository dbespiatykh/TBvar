from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider

NCBI = NCBIRemoteProvider(email=config["NCBI"]["email"])


rule download_genome:
    input:
        NCBI.remote("NC_000962.3.fasta", db="nuccore"),
    output:
        "ref/NC_000962.3.fa",
    log:
        "logs/ncbi/ref_dwn.log",
    run:
        shell("cat {input} > {output}")


rule samtools_genome_index:
    input:
        "ref/NC_000962.3.fa",
    output:
        "ref/NC_000962.3.fa.fai",
    log:
        "logs/samtools/ref_index.log",
    wrapper:
        "v1.7.1/bio/samtools/faidx"


rule create_dict:
    input:
        "ref/NC_000962.3.fa",
    output:
        "ref/NC_000962.3.dict",
    log:
        "logs/picard/create_dict.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.7.1/bio/picard/createsequencedictionary"


rule bwa_index:
    input:
        "ref/NC_000962.3.fa",
    output:
        idx=multiext("ref/NC_000962.3.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa/ref_index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.7.1/bio/bwa/index"

rule barcode_snps:
    input:
        "results/all.snps.pass.txt",
        config["files"]["barcodes"],
    output:
        "results/samples.barcoded.snps.tsv",
    log:
        "logs/barcoding/barcoding_snps.log",
    conda:
        "../envs/barcoding.yaml"
    script:
        "../scripts/barcoding_snps.py"

rule barcode_levels:
    input:
        "results/all.snps.pass.txt",
        config["files"]["levels"],
    output:
        "results/samples.barcoded.levels.tsv",
    log:
        "logs/barcoding/barcoding_levels.log",
    conda:
        "../envs/barcoding.yaml"
    script:
        "../scripts/barcoding_levels.py"
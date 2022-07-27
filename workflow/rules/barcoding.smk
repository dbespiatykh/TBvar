## Make vcf stats table
rule barcode_samples:
    input:
        "results/all.snps.pass.txt",
        config["files"]["barcodes"],
    output:
        "results/all.snps.pass.barcoded.tsv",
    log:
        "logs/barcoding/barcoding.log",
    conda:
        "../envs/barcoding.yaml"
    script:
        "../scripts/barcoding.py"

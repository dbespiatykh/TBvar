rule barcode_snps:
    input:
        "results/all.snps.pass.txt",
        config["files"]["barcodes"],
    output:
        report(
            "results/samples.barcoded.snps.tsv",
            caption="../report/barcoded_snps.rst",
            category="Barcoding",
        ),
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
        report(
            "results/samples.barcoded.levels.tsv",
            caption="../report/barcoded_levels.rst",
            category="Barcoding",
        ),
    log:
        "logs/barcoding/barcoding_levels.log",
    conda:
        "../envs/barcoding.yaml"
    script:
        "../scripts/barcoding_levels.py"

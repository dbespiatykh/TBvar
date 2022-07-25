## Make vcf stats table
rule barcode_samples:
    input:
        "Results/all.snps.pass.txt",
        "barcodes/barcodes.tsv"
    output:
        "Results/all.snps.pass.barcoded.tsv",
    conda:
        "../envs/barcoding.yaml"
    script:
        "../scripts/barcoding.py"
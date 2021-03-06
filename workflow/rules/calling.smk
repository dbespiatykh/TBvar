## Call variants
rule gatk_haplotype_caller:
    input:
        bam="BAM/{sample}.bam",
        idx="BAM/{sample}.bam.bai",
        ref="ref/NC_000962.3.fa",
        fai="ref/NC_000962.3.fa.fai",
        dic="ref/NC_000962.3.dict",
    output:
        gvcf="VCF/gVCF/{sample}.g.vcf.gz",
    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    params:
        extra="-ploidy 1 -mbq 20",
    threads: config["GATK"]["hapcall"]["threads"]
    resources:
        mem_mb=config["GATK"]["hapcall"]["memory"],
    wrapper:
        "v1.7.1/bio/gatk/haplotypecaller"


## Import gVCFs to GenomicsDB
rule gatk_genomics_db_import:
    input:
        gvcfs=expand("VCF/gVCF/{sample}.g.vcf.gz", sample=samples.index),
    output:
        db=temp(directory("VCF/db")),
    log:
        "logs/gatk/genomicsdbimport.log",
    params:
        intervals="NC_000962.3",
        db_action="create",
        extra="--batch-size 100",
    resources:
        mem_mb=config["GATK"]["gendb"]["memory"],
    wrapper:
        "v1.7.1/bio/gatk/genomicsdbimport"


## Genotype Variants
rule gatk_genotype_gvcfs:
    input:
        genomicsdb="VCF/db",
        ref="ref/NC_000962.3.fa",
    output:
        vcf="VCF/all.vcf.gz",
    log:
        "logs/gatk/genotypegvcfs.log",
    params:
        extra="-ploidy 1",
    resources:
        mem_mb=config["GATK"]["genotype"]["memory"],
    wrapper:
        "v1.7.1/bio/gatk/genotypegvcfs"


## Filter SNPs
rule gatk_filter_variants:
    input:
        vcf="VCF/all.vcf.gz",
        ref="ref/NC_000962.3.fa",
    output:
        vcf="VCF/all.filtered.vcf.gz",
    log:
        "logs/gatk/filter/snps.filter.log",
    params:
        filters={"Failfilter": "QD < 2.0 || DP < 10 || FS > 60.0 || MQ < 40.0"},
    resources:
        mem_mb=config["GATK"]["filter"]["memory"],
    wrapper:
        "v1.7.1/bio/gatk/variantfiltration"


## Select only PASS SNPs
rule gatk_select_variants:
    input:
        vcf="VCF/all.filtered.vcf.gz",
        ref="ref/NC_000962.3.fa",
    output:
        vcf="VCF/all.snps.pass.vcf.gz",
    log:
        "logs/gatk/select/snps.pass.log",
    params:
        extra="--exclude-filtered --select-type-to-include SNP",
    resources:
        mem_mb=config["GATK"]["select"]["memory"],
    wrapper:
        "v1.7.1/bio/gatk/selectvariants"


## Refactor Variants to Table
rule gatk_variants_to_table:
    input:
        vcf="VCF/all.snps.pass.vcf.gz",
    output:
        tab="results/all.snps.pass.txt",
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/gatk/vartotable/vars2table.log",
    shell:
        "gatk VariantsToTable \
        -V {input.vcf} \
        -O {output.tab} \
        -F POS \
        -F REF \
        -GF GT 2> {log}"

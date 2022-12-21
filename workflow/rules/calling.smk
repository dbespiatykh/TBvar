## Call variants
rule gatk_haplotype_caller:
    input:
        bam="results/BAM/{sample}.bam",
        idx="results/BAM/{sample}.bam.bai",
        ref="resources/ref/NC_000962.3.fa",
        fai="resources/ref/NC_000962.3.fa.fai",
        dic="resources/ref/NC_000962.3.dict",
    output:
        gvcf="results/VCF/gVCF/{sample}.g.vcf.gz",
    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    params:
        extra="-ploidy 1 -mbq 20",
    threads: config["GATK"]["hapcall"]["threads"]
    resources:
        mem_mb=config["GATK"]["hapcall"]["memory"],
    wrapper:
        "v1.21.0/bio/gatk/haplotypecaller"


## Import gVCFs to GenomicsDB
rule gatk_genomics_db_import:
    input:
        gvcfs=[
            expand("results/VCF/gVCF/{sample}.g.vcf.gz", sample=samples.index),
            config["files"]["dummy"],
        ],
    output:
        db=temp(directory("results/VCF/db")),
    log:
        "logs/gatk/genomicsdbimport.log",
    params:
        intervals="NC_000962.3",
        db_action="create",
        extra="--batch-size 100",
    resources:
        mem_mb=config["GATK"]["gendb"]["memory"],
    wrapper:
        "v1.21.0/bio/gatk/genomicsdbimport"


## Genotype Variants
rule gatk_genotype_gvcfs:
    input:
        genomicsdb="results/VCF/db",
        ref="resources/ref/NC_000962.3.fa",
        fai="resources/ref/NC_000962.3.fa.fai",
        dic="resources/ref/NC_000962.3.dict",
    output:
        vcf="results/VCF/all.vcf.gz",
    log:
        "logs/gatk/genotypegvcfs.log",
    params:
        extra="-ploidy 1",
    resources:
        mem_mb=config["GATK"]["genotype"]["memory"],
    wrapper:
        "v1.21.0/bio/gatk/genotypegvcfs"


## Filter SNPs
rule gatk_filter_variants:
    input:
        vcf="results/VCF/all.vcf.gz",
        ref="resources/ref/NC_000962.3.fa",
        fai="resources/ref/NC_000962.3.fa.fai",
        dic="resources/ref/NC_000962.3.dict",
    output:
        vcf="results/VCF/all.filtered.vcf.gz",
    log:
        "logs/gatk/filter/snps.filter.log",
    params:
        filters={"Failfilter": "QD < 2.0 || DP < 10 || FS > 60.0 || MQ < 40.0"},
    resources:
        mem_mb=config["GATK"]["filter"]["memory"],
    wrapper:
        "v1.21.0/bio/gatk/variantfiltration"


## Split deletions and snps from single position
rule gatk_left_align_and_trim_variants:
    input:
        vcf="results/VCF/all.filtered.vcf.gz",
        ref="resources/ref/NC_000962.3.fa",
        fai="resources/ref/NC_000962.3.fa.fai",
        dic="resources/ref/NC_000962.3.dict",
    output:
        vcf="results/VCF/all.filtered.trimmed.vcf.gz",
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/gatk/leftalingandtrim/leftalingandtrim.log",
    shell:
        "gatk LeftAlignAndTrimVariants \
        -V {input.vcf} \
        -O {output.vcf} \
        -R {input.ref} \
        --split-multi-allelics 2> {log}"


## Select only PASS SNPs
rule gatk_select_variants:
    input:
        vcf="results/VCF/all.filtered.trimmed.vcf.gz",
        ref="resources/ref/NC_000962.3.fa",
        fai="resources/ref/NC_000962.3.fa.fai",
        dic="resources/ref/NC_000962.3.dict",
    output:
        vcf=report(
            "results/VCF/all.snps.pass.vcf.gz",
            caption="../report/vcf.rst",
            category="SNPs",
        ),
    log:
        "logs/gatk/select/snps.pass.log",
    params:
        extra="--exclude-filtered --select-type-to-include SNP --exclude-sample-name Dummy",
    resources:
        mem_mb=config["GATK"]["select"]["memory"],
    wrapper:
        "v1.21.0/bio/gatk/selectvariants"


## Refactor Variants to Table
rule gatk_variants_to_table:
    input:
        vcf="results/VCF/all.snps.pass.vcf.gz",
        intervals=config["files"]["intervals"],
    output:
        tab=report(
            "results/all.snps.pass.txt",
            caption="../report/tab.rst",
            category="SNPs",
        ),
    conda:
        "../envs/gatk4.yaml"
    log:
        "logs/gatk/vartotable/vars2table.log",
    shell:
        "gatk VariantsToTable \
        -V {input.vcf} \
        -O {output.tab} \
        --intervals {input.intervals} \
        --interval-padding 1 \
        -F POS \
        -F REF \
        -GF GT 2> {log}"

Reads were mapped to the referrence *M. tuberculosis* H37Rv genome (RefSeq: `NC_000962.3`_) using `BWA MEM`_ then
mapped reads were sorted with `SAMtools`_ and duplicate reads were removed using `Picard MarkDuplicates`_.
Variants in `gVCF`_ format were called with `GATK HaplotypeCaller`_ and imported into GenomicsDB before joint genotyping using
`GATK GenomicsDBImport`_. Then joint genotyping was performed with `GATK GenotypeGVCFs`_ tool. Called variants were filtered 
with `GATK VariantFiltration`_ using the following filter expression: `QD`_ < 2.0 || `DP`_ < 10 || `FS`_ > 60.0 || `MQ`_ < 40.0. 
Multiallics variants were then split into biallelics, left aligned and trimmed with `GATK LeftAlignAndTrimVariants`_. Only SNPs that 
passed filtering were subsetted using `GATK SelectVariants`_ and subsequently transformed into tab-delimited format with `GATK VariantsToTable`_.


At the final step lineage barcoding for each sample was done using `Python`_ scripts.


.. _NC_000962.3: https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3/
.. _BWA MEM: https://github.com/lh3/bwa
.. _SAMtools: https://samtools.sourceforge.net/
.. _Picard MarkDuplicates: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
.. _gVCF: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
.. _GATK HaplotypeCaller: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
.. _GATK GenomicsDBImport: https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport
.. _GATK GenotypeGVCFs: https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs
.. _GATK VariantFiltration: https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration#--filter-expression
.. _GATK SelectVariants: https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
.. _QD: https://gatk.broadinstitute.org/hc/en-us/articles/360036818051-QualByDepth
.. _DP: https://gatk.broadinstitute.org/hc/en-us/articles/4404604815643-Coverage
.. _FS: https://gatk.broadinstitute.org/hc/en-us/articles/360040096152-FisherStrand
.. _MQ: https://gatk.broadinstitute.org/hc/en-us/articles/360036714651-RMSMappingQuality
.. _GATK LeftAlignAndTrimVariants: https://gatk.broadinstitute.org/hc/en-us/articles/360037225872-LeftAlignAndTrimVariants
.. _GATK VariantsToTable: https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable
.. _Python: https://www.python.org/

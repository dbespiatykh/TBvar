import pandas as pd
import pandas.io.formats.excel

pandas.io.formats.excel.ExcelFormatter.header_style = None

vcf_report = pd.read_csv(
    snakemake.input[0], sep='\t', header=0, names=['Sample', 'SNPs']
).sort_values(by=['Sample'])
fastqc_report = pd.read_csv(
    snakemake.input[1],
    sep='\t',
    usecols=[
        "Sample",
        "Total Sequences",
        "Sequence length",
        "%GC",
        "avg_sequence_length",
    ],
).sort_values(by=['Sample'])
samtools_report = pd.read_csv(
    snakemake.input[2],
    sep='\t',
    usecols=[
        "Sample",
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "reads_paired",
        "reads_mapped_percent",
        "reads_unmapped_percent",
        "average_quality",
    ],
).sort_values(by=['Sample'])
barcodes = pd.read_csv(snakemake.input[3], sep='\t').sort_values(by=['sample'])
mosdepth = pd.read_csv(
    snakemake.input[4], sep='\t', header=0, names=['Sample', 'Mean depth']
).sort_values(by=['Sample'])
spoligotypes = pd.read_csv(snakemake.input[5], sep='\t').sort_values(by=['StrainID'])
tbprofiler = pd.read_csv(snakemake.input[6], sep='\t').sort_values(by=['sample'])


def pe_se(samtools_report):
    if samtools_report['reads_paired'] > 0:
        return "paired_end"
    elif samtools_report['reads_paired'] == 0:
        return "single_end"


samtools_report['library_type'] = samtools_report.apply(pe_se, axis=1)
vcf_fastq = pd.merge(vcf_report, fastqc_report, on='Sample', how="left")
all_stats = pd.merge(vcf_fastq, samtools_report, on='Sample', how="left")
fin_stats = pd.merge(all_stats, mosdepth, on='Sample', how="left")

writer = pd.ExcelWriter(snakemake.output[0], engine='xlsxwriter')

fin_stats.to_excel(writer, sheet_name='Stats', index=False)
barcodes.to_excel(writer, sheet_name='Barcodes', index=False)
spoligotypes.to_excel(writer, sheet_name='Spoligotypes', index=False)
tbprofiler.to_excel(writer, sheet_name='Drug resistance (TBProfiler)', index=False)

writer.save()


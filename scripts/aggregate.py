import os
import pandas as pd

vcf_report = pd.read_csv(snakemake.input[0], sep='\t').sort_values(by=['SampleName'])
fastqc_report = pd.read_csv(snakemake.input[1], sep='\t').sort_values(by=['Sample'])
samtools_report = pd.read_csv(snakemake.input[2], sep='\t').sort_values(by=['Sample'])
picard_report = pd.read_csv(snakemake.input[3], sep='\t').sort_values(by=['Sample'])
barcodes = pd.read_csv(snakemake.input[4], sep='\t').sort_values(by=['sample'])
spoligotypes = pd.read_csv(snakemake.input[5], sep='\t').sort_values(by=['StrainID'])
tbprofiler = pd.read_csv(snakemake.input[6], sep='\t').sort_values(by=['sample'])

writer = pd.ExcelWriter(snakemake.output[0], engine='xlsxwriter')

vcf_report.to_excel(writer, sheet_name='Variants stats', index=False)
fastqc_report.to_excel(writer, sheet_name='Reads stats', index=False)
samtools_report.to_excel(writer, sheet_name='Mapping stats', index=False)
picard_report.to_excel(writer, sheet_name='Duplication stats', index=False)
barcodes.to_excel(writer, sheet_name='Barcodes', index=False)
spoligotypes.to_excel(writer, sheet_name='Spoligotypes', index=False)
tbprofiler.to_excel(writer, sheet_name='Drug resistance (TBProfiler)', index=False)

writer.save()

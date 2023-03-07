# General settings

Following the explanations in the `config.yml`, you should modify it according to your needs.

# Samples sheet

The location of this sheet must be specified in the `config.yml`.

It should be formatted like this:

| Run_accession | R1                             | R2                             |
| ------------- | ------------------------------ | ------------------------------ |
| SRR2024996    | /path/to/SRR2024996_1.fastq.gz | /path/to/SRR2024996_2.fastq.gz |
| SRR2024925    | /path/to/SRR2024925_1.fastq.gz | /path/to/SRR2024925_2.fastq.gz |
| SRR12882189   | /path/to/SRR12882189.fastq.gz  |                                |

**Run_accession** - Run accession number or sample name;\
**R1** - Path to the first read pair;\
**R2** - Path to the second read pair.

If samples are PE:

- Both `R1` and `R2`should be specified with reads paths.

If samples are SE:

- Only `R1` column should have the path to the SE reads file.

You can use the following code to create `samples.tsv` file:

```bash
for f in path/to/fastq/dir/*_1.fastq.gz; do printf '%s\t%s\t%s\n' "$(basename "$f" _1.fastq.gz)" "$f" "${f%_1.fastq.gz}_2.fastq.gz"; done | sed -E -e '1iRun_accession\tR1\tR2' > samples.tsv
```

Replace `path/to/fast/dir/` with the path to your directory containing FASTQ files, and `_1` and `_2` if you read suffixes are different.

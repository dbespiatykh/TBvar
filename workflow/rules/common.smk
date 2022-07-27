import pandas as pd


configfile: "config/config.yml"


samples = (
    pd.read_csv(
        config["files"]["samples"],
        sep="\t",
        dtype=str,
    )
    .set_index(["Run_accession"], drop=False)
    .sort_index()
)


def get_fastq(wildcards):
    sample = samples.loc[wildcards.sample]

    if pd.isna(sample["R2"]):
        return [sample["R1"]]
    else:
        return [sample["R1"], sample["R2"]]

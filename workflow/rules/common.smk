import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("7.18.2")


configfile: "config/config.yml"


report: "../report/workflow.rst"


validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(
        config["files"]["samples"],
        sep="\t",
        dtype=str,
    )
    .set_index(["Run_accession"], drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")


def get_fastq(wildcards):
    sample = samples.loc[wildcards.sample]

    if pd.isna(sample["R2"]):
        return [sample["R1"]]
    else:
        return [sample["R1"], sample["R2"]]

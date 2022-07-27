import pandas as pd

configfile: "config/config.yml"
samples = pd.read_csv(config["smp_file"], dtype=str, sep="\t").set_index(
    "Run_accession", drop=False
)

sample = samples.index.to_list()[0]


def get_fastq(wildcard):
    return samples.loc[wildcard, ["R1", "R2"]].dropna().to_dict()

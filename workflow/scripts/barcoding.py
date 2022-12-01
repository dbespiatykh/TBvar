#!/usr/bin/env python3

import numpy as np
import pandas as pd

chunks = pd.read_csv(snakemake.input[0], sep="\t", chunksize=1000)

barcodes = pd.read_csv(snakemake.input[1], sep="\t")

chunk_list = []

for chunk in chunks:
    result = chunk.melt(
        chunk[chunk.columns[chunk.columns.isin(["POS", "REF"])]],
        var_name="Sample",
        value_name="ALT",
    )
    result["Sample"] = result["Sample"].str.replace(r".GT", r"", regex=False)
    result = result.drop(result[result["REF"] == result["ALT"]].index)
    chunk_list.append(result)

df = pd.concat(chunk_list).reset_index(drop=True)

pos_b = barcodes["POS"].values
ref_b = barcodes["REF"].values
alt_b = barcodes["ALT"].values
pos_e = df["POS"].values[:, None]
ref_e = df["REF"].values[:, None]
alt_e = df["ALT"].values[:, None]

df["lineage"] = np.dot(
    np.logical_and.reduce(
        [
        np.equal(pos_b, pos_e),
        np.equal(ref_b, ref_e),
        np.equal(alt_b, alt_e)
        ]
    ),
    barcodes["lineage"],
)

df = df.drop(["REF", "ALT", "POS"], axis=1)

df = (
    df.groupby(["Sample"])
    .agg(lambda x: ",".join(set(x)))
    .reset_index()
    .reindex(columns=df.columns)
)

df = df.join(df["lineage"].str.get_dummies(sep=",").astype(int))

df = df.drop(["lineage"], axis=1).replace({1: "+", 0: np.nan})

cols = barcodes.iloc[:, 0].to_list()

cols.insert(0, "Sample")

df.columns = pd.CategoricalIndex(df.columns, categories=cols, ordered=True)

df = df.sort_index(axis=1)

df.to_csv(snakemake.output[0], sep="\t", index=False)

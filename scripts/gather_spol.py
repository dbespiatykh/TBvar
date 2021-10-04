import os
import pandas as pd

df = pd.concat(
    [
        pd.read_csv(
            fp, sep="\t", names=["path", "bin", "oct"], dtype={"bin": str, "oct": str}
        ).assign(Sample=os.path.basename(fp).split(".SpoTyping")[0])
        for fp in snakemake.input
    ]
)
df.drop(columns=["path", "oct"], inplace=True)
df = df[["Sample", "bin"]]
df.to_csv(snakemake.output[0], sep="\t", index=False, header=False)

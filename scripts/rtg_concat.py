# %%
import os
import pandas as pd

# %%
with open(snakemake.input[0], 'r') as file:
    data = file.read().split('\n\n')
data = [i.split('\n') for i in data]

# %%
df = pd.DataFrame(data)

df=df[[1,2]]

df.rename(columns={1: "SampleName", 2: "SNPs"}, inplace=True)

# %%
df = df[['SampleName', 'SNPs']]

# %%
df.to_csv(snakemake.output[0], sep='\t', index=False)

# %%

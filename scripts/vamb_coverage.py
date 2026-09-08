"""Reformat the clean-coverage table into Vamb's expected abundance TSV."""

import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t")
df = df.rename(columns={"contigName": "contigname"})
df = df.loc[:, (df != 0).any(axis=0)]
df.to_csv(snakemake.output[0], sep="\t", index=False, header=True)

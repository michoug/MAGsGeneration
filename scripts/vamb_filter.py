"""Restrict the taxonomy table to contigs that survived the Vamb length filter."""

import pandas as pd

ids = []
with open(snakemake.input.fasta) as file:
    for line in file:
        if line.startswith(">"):
            ids.append(line[1:].strip())

df = pd.read_csv(snakemake.input.tax, sep="\t")
filtered_df = df[df.iloc[:, 0].isin(ids)]
filtered_df.to_csv(snakemake.output[0], sep="\t", index=False, header=True)

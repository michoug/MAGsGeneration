"""Fill in "unknown" taxonomy for contigs with no mmseqs prediction."""

import pandas as pd

contigs = []
with open(snakemake.input.fasta) as f:
    for line in f:
        if line.startswith(">"):
            contigs.append(line[1:].strip().split()[0])

tax_df = pd.read_csv(
    snakemake.input.tax, sep="\t", header=None, names=["contigs", "predictions"]
)
all_df = pd.DataFrame({"contigs": contigs})
merged = all_df.merge(tax_df, on="contigs", how="left")
merged["predictions"] = merged["predictions"].fillna("unknown")
merged.to_csv(snakemake.output[0], sep="\t", index=False, header=True)

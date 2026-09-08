"""Compute the percentage of eukaryotic vs prokaryotic contigs from DeepMicroClass predictions."""

import os

import pandas as pd

df = pd.read_csv(snakemake.input.pred, sep="\t", header=0)
choice_labels = {
    1: "Eukaryote",
    2: "Eukaryote",
    3: "Prokaryote",
    4: "Prokaryote",
    5: "Prokaryote",
}
df["pred"] = df["best_choice"].map(choice_labels)
count = pd.DataFrame(df.groupby("pred")["pred"].count())
count["percent"] = (count["pred"] / count["pred"].sum()) * 100
count_trans = count.transpose()
count_trans["bin"] = os.path.basename(snakemake.input.pred)
count_trans["bin"] = count_trans["bin"].str.replace("_pred_one-hot_hybrid.tsv", "")
count_trans_edited = count_trans.drop(labels="pred")
count_trans_edited.to_csv(snakemake.output.percent, sep="\t", header=True)

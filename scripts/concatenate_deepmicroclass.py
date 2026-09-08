"""Collate per-bin DeepMicroClass percentages and list bins predicted as eukaryotic."""

import contextlib

import pandas as pd

with (
    open(snakemake.log.out, "w") as f_log, contextlib.redirect_stdout(f_log),
    contextlib.redirect_stderr(f_log),
):
    dfs = []
    for file in snakemake.input.all_files:
        df = pd.read_csv(file, sep="\t", header=0)
        print(f"Loaded file: {file}, shape: {df.shape}")
        dfs.append(df)
    concatenated = pd.concat(dfs, ignore_index=True)
    print(f"Concatenated DataFrame shape: {concatenated.shape}")
    concatenated.to_csv(snakemake.output.collated, sep="\t", index=False)
    print(f"Collated output written to: {snakemake.output.collated}")
    eukaryotic_bins = concatenated[concatenated["Eukaryote"] >= 80]
    eukaryotic_bins[["bin"]].to_csv(
        snakemake.output.final_list, sep="\t", index=False, header=False
    )
    print(f"Eukaryotic bins written to: {snakemake.output.final_list}")

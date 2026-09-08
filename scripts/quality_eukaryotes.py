"""Extract quality metrics for eukaryotic MAGs from the BUSCO MultiQC report."""

from pathlib import Path

import pandas as pd

df = pd.read_csv(snakemake.input.stats, sep="\t")
src = pd.read_csv(snakemake.input.source, sep="\t")
required_cols = {
    "Sample",
    "busco-complete",
    "busco-complete_single_copy",
    "busco-fragmented",
    "busco-missing",
}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    raise ValueError(f"Missing required columns: {sorted(missing_cols)}")
required_src_cols = {"Sample Name", "Source"}
missing_src_cols = required_src_cols - set(src.columns)
if missing_src_cols:
    raise ValueError(
        f"Missing required columns in sources file: {sorted(missing_src_cols)}"
    )
busco_total = (
    df["busco-complete_single_copy"]
    + df["busco-fragmented"]
    + df["busco-missing"]
)
if "busco-complete_duplicated" in df.columns:
    duplicated = df["busco-complete_duplicated"]
else:
    duplicated = (df["busco-complete"] - df["busco-complete_single_copy"]).clip(
        lower=0
    )
src["sample_key"] = src["Sample Name"].str.strip()
src["genome"] = (
    src["Source"]
    .str.extract(r"/eukaryotic_mags_busco/([^/]+)/", expand=False)
    .fillna(src["Source"].str.extract(r"/BUSCO_([^/]+)\.fa(?:/|$)", expand=False))
)
sample_to_genome = (
    src.dropna(subset=["genome"])
    .drop_duplicates(subset=["sample_key"])
    .set_index("sample_key")["genome"]
)
df["sample_key"] = df["Sample"].str.strip()
df["genome"] = df["sample_key"].map(sample_to_genome)
fallback = (
    df["sample_key"]
    .str.replace(r"^short_summary\.specific\.eukaryota_odb12\.", "", regex=True)
    .str.replace(r"^BUSCO_", "", regex=True)
    .str.replace(r"\.fa(?:\.txt)?$", "", regex=True)
)
df["genome"] = df["genome"].fillna(fallback)
df["completeness"] = (
    (df["busco-complete"] + df["busco-fragmented"]) * 100 / busco_total
)
df["contamination"] = duplicated * 100 / busco_total
result = (
    df[["genome", "completeness", "contamination"]]
    .sort_values("completeness", ascending=False)
    .drop_duplicates(subset="genome", keep="first")
)
mags_folder = Path(snakemake.input.mags_done).parent
all_genomes = {p.stem for p in mags_folder.glob("*.fa")}
present_genomes = set(result["genome"].astype(str))
missing_genomes = sorted(all_genomes - present_genomes)
if missing_genomes:
    missing_df = pd.DataFrame(
        {
            "genome": missing_genomes,
            "completeness": 0.0,
            "contamination": 0.0,
        }
    )
    result = pd.concat([result, missing_df], ignore_index=True)
result = result.sort_values("genome")
result.to_csv(snakemake.output[0], index=False)

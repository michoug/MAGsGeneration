# MAGs Generation Pipeline

## Pipeline description

There are two main options of this pipeline, either you have assemblies and want to bin them, or you have raw reads and want to go from reads to MAGs. 

### Mandatory steps

- The reads are trimmed with [trim_galore](https://github.com/FelixKrueger/TrimGalore) and human reads are removed with [strobealign](https://github.com/ksahlin/strobealign) and [samtools](https://github.com/samtools/samtools) across the human genome (hg38, GCF_000001405.40). Duplicate reads are then removed with [fastp](https://github.com/OpenGene/fastp).
- The trimmed reads are then assembled with [megahit](https://github.com/voutcn/megahit) (`runAll` pipeline) or not (`alreadyAssemble` pipeline). Contigs shorter than 1.5 kb are filtered out with [seqkit](https://github.com/shenwei356/seqkit).
- Contig taxonomy is assigned using [MMseqs2](https://github.com/soedinglab/MMseqs2) against the GTDB database. This taxonomy is used as input for the taxonomy-aware binning step.
- The coverage is calculated with [strobealign](https://github.com/ksahlin/strobealign) by mapping all reads across all assemblies.
- Bins are generated with [Comebin](https://github.com/ziyewang/COMEBin), [Concoct](https://github.com/BinPro/CONCOCT), [MetaCAT](https://github.com/liu-congcong/MetaCAT), [Metabat2](https://bitbucket.org/berkeleylab/metabat), [Rosella](https://github.com/rhysnewell/rosella), [Semibin2](https://github.com/BigDataBiology/SemiBin) and [Taxvamb](https://github.com/RasmussenLab/vamb) per sample.
- Each binner's raw bins are scored independently with [CheckM2](https://github.com/chklovski/CheckM2) and ranked using the methodology of [Han, Wang & Zhu (Nat Commun 16:2865, 2025)](https://doi.org/10.1038/s41467-025-57957-6): binners are grouped into three technology families (`classical` feature-based clustering - Metabat2, Concoct; `ml_clustering` - Rosella, MetaCAT; `deep_representation` - Comebin, Semibin2, Taxvamb), scored per sample on the count of bins reaching medium/near-complete/strict-near-complete quality tiers, and the best-scoring binner from each family (up to `binners.selection.n`, default 3) is kept. This mirrors the paper's finding that refinement tools perform best when fed a handful of complementary binners rather than all of them; it approximates the paper's HQ tier (which also requires rRNA/tRNA gene presence) with a stricter completeness/contamination tier to avoid adding barrnap/tRNAscan-SE to the pipeline. Set `binners.selection.enabled: false` in the config to skip scoring and pass all seven binners through instead.
- [DAS Tool](https://github.com/cmks/DAS_Tool) is then used to select the best bins from the selected binners per sample.
- The DAS Tool bins are refined using [Rosella](https://github.com/rhysnewell/rosella).
- All bins from all samples are then combined in a single folder and [CheckM2](https://github.com/chklovski/CheckM2) is used to estimate quality (completeness and contamination).
- All the bins that pass a quality threshold of 70% completeness and <10% contamination are dereplicated with [Galah](https://github.com/wwood/galah).
- [CheckM2](https://github.com/chklovski/CheckM2) and [GtdbTk](https://github.com/Ecogenomics/GTDBTk) are then used to estimate the quality and taxonomy of the final MAGs.
- Read mapping of all samples against the final sets of MAGs are done using [CoverM](https://github.com/wwood/CoverM).

### Optional steps

- Annotations of the MAGs is done with [pyrodigal](https://github.com/althonos/pyrodigal) and [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper/)
- [MGThermometer](https://doi.org/10.1101/2022.07.14.499854) is used to measure the `optimal growth rate` based on the relative abundance of `FIVYWREL` aminoacids
  - Optimal growth rate is measured as follows,
```math
OGT = 937 * F_{IVYWREL} − 335
```
- [SynTracker](https://github.com/leylabmpi/SynTracker) is used to track strains across samples based on synteny
- [SingleM](https://github.com/wwood/singlem) is a read-based software used to profile the samples based on single copy marker genes (mostly ribosomal) and assess how much of the diversity is captured in the MAGs
- Putative eukaryotic bins are first re-assembled and refined with [remag](https://github.com/deng-lab/remag), then filtered based on size (>2.5 Mb) and classified with [DeepMicroClass](https://github.com/chengsly/DeepMicroClass), which distinguishes prokaryotic from eukaryotic contigs. Bins with >80% of eukaryotic contigs are considered eukaryotic bins. Their completeness and contamination are then estimated with [Busco](https://busco.ezlab.org/) and dereplicated with [Galah](https://github.com/wwood/galah).

## Setup

### Pixi

This pipeline uses **[Pixi](https://prefix.dev/pixi)** to manage all tool environments. Pixi automatically installs isolated per-tool environments from `pixi.toml` — no manual conda setup required.

```bash
# install pixi
curl -fsSL https://pixi.sh/install.sh | sh
```

Getting the repository
```bash
git clone https://github.com/michoug/MAGsGeneration.git
cd MAGsGeneration
```

Install the default Snakemake environment

```bash
pixi install
```

All per-step tool environments are installed automatically when the pipeline first runs them, or you can pre-install everything at once:

```bash
pixi install -a --run-post-link-scripts
```

### Run Setup

* Place your reads (e.g. `sample_R1_001.fastq.gz` and `sample_R2_001.fastq.gz` files) in a `reads` folder
* Eventually, place your assemblies (e.g. `sample.fna` files) in an `assemblies` folder
* Fill the `config/sample_list.txt` file with the name of your samples (e.g. `{sample}_R2_001.fastq.gz`)
* Modify the `config/config.yaml` file to change the different paths and eventually the different options, especially `work_dir`, `results_dir`, `reads_dir` and `assemblies_dir`

## Running the pipeline

### Via Python script

This Python script provides a flexible launcher for running the NOMIS MAGs analysis pipelines with Snakemake. It supports **local execution**, **Slurm cluster execution**, **optional steps**, **continuation of failed runs**, and **dry-run previews**.
The pipeline works with **[Pixi](https://prefix.dev/pixi)** to manage environments for each step. It installs **Pixi** if it's missing and then installs all required environments from `pixi.toml`.

#### Usage

##### Local Execution

```bash
pixi r python MAGsGeneration.py -p <pipeline> -c <config.yaml> --output <results_path> [options]
```

##### Slurm Execution

```bash
pixi r python MAGsGeneration.py -p <pipeline> -c <config.yaml> --output <results_path> --slurm [options]
```

##### Dry-Run

Preview the steps without executing:

```bash
pixi r python MAGsGeneration.py -c <results_path>/config_run.yaml --continue --dryrun --cores 4
```

##### Continue Failed Run

Resume unfinished or failed steps:

```bash
pixi r python MAGsGeneration.py -c <results_path>/config_run.yaml --continue --cores 4
```

---

#### Arguments

| Argument         | Description                                                                   |
| ---------------- | ----------------------------------------------------------------------------- |
| `-p, --pipeline` | Pipeline to run (`runAll`, `alreadyAssemble`). Not required for continuation. |
| `-c, --config`   | Path to the base or previous run YAML config.                                 |
| `--steps`        | Optional steps to include (defaults to all optional steps).                   |
| `--output`       | Path to the results folder (full or relative path).                           |
| `--dryrun`       | Perform a dry-run to preview steps.                                           |
| `--profile`      | Base Snakemake profile directory (default: `profile/`).                       |
| `--slurm`        | Use the Slurm profile for execution.                                          |
| `--use-conda`    | Enable `--use-conda` in Snakemake.                                            |
| `--conda-prefix` | Conda environment prefix path.                                                |
| `--cores`        | Number of cores to use.                                                       |
| `--continue`     | Continue previous run using the same results folder.                          |

---

#### Optional Steps

By default, all steps are included. The following **optional steps** can be toggled via `--steps`:

* `eggnog`
* `eukaryotes_filter`
* `mgthermometer`
* `singlem_all`
* `syntracker`

Mandatory steps are always executed. The launcher prints a summary before running:

```
🚀 Steps to be executed:
  Mandatory steps: ['install_software', 'trim_reads', 'assembly_all', 'checkm']
  Optional steps:  ['syntracker', 'eggnog']
```

---

#### Output and Metadata

* The **results folder** is specified via `--output`.

* The **run configuration** is saved as `config_run.yaml` inside the results folder.

* A **metadata file** `run_info.txt` is generated including:

  * Timestamp
  * Snakefile used
  * Profile path
  * Steps executed
  * Dry-run status

* Using the same `--output` folder allows **continuation** of failed runs or dry-runs without creating new directories.


### Via Snakemake

You can also invoke Snakemake directly without the Python wrapper. In this case you are responsible for preparing the config file and managing the working directory yourself.

#### Local execution (without Slurm)

```bash
pixi r snakemake -s workflow/runAll.smk --configfile config/config.yaml --cores 4
```

#### Dry run (preview steps without executing)

```bash
pixi r snakemake -s workflow/runAll.smk --configfile config/config.yaml --dry-run
```

#### With Slurm

```bash
pixi r snakemake -s workflow/runAll.smk --configfile config/config.yaml --profile profile/
```

#### On already-assembled data

```bash
pixi r snakemake -s workflow/alreadyAssemble.smk --configfile config/config.yaml --cores 4
```

#### Unlock the working directory (if a previous run was interrupted)

```bash
pixi r snakemake -s workflow/runAll.smk --configfile config/config.yaml --unlock
```

#### Resume an incomplete run

```bash
pixi r snakemake -s workflow/runAll.smk --configfile config/config.yaml --cores 4 --rerun-incomplete
```


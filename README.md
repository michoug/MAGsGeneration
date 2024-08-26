# MAGs Generation Pipeline

## Pipeline description

- This pipeline expect contigs samples (`fa` files) in the `assembly_dir` and their respective reads (`_1_val_1.fq.gz` or `_2_val_2.fq.gz` files) in the `reads_dir`.
- The coverage is calculated with [fairy](https://github.com/bluenote-1577/fairy) for binning that speeds up the process greatly
- Bins are generated with [Metabat2](https://bitbucket.org/berkeleylab/metabat), [Metabinner](https://github.com/ziyewang/MetaBinner) and [Rosella](https://github.com/rhysnewell/rosella) per sample.
- [DAS Tool](https://github.com/cmks/DAS_Tool) is then used to generate a non redudant set of bins per samples.
- These redundant sets are refined using [Rosella](https://github.com/rhysnewell/rosella).
- All bins from all samples are then combined in a single folder and [Checkm2](https://github.com/chklovski/CheckM2) is used to estimate quality (completeness and contamination)
- All the bins that passes a quality threshold of 70% completeness and <10% contamination are dereplicated with [Galah](https://github.com/wwood/galah).
- [CheckM2](https://github.com/chklovski/CheckM2) and [GtdbTk](https://github.com/Ecogenomics/GTDBTk) are then used to estimate the quality and taxonomy of the final MAGs.
- Read mapping of all samples against the final sets of MAGs are done using [CoverM](https://github.com/wwood/CoverM). 
- Annotations of the MAGs is done with [pyrodigal](https://github.com/althonos/pyrodigal) and [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper/)
- [MGThermometer](https://doi.org/10.1101/2022.07.14.499854) is used to measure the `optimal growth rate` based on the relative abundance of `FIVYWREL` aminoacids
  - Optimal growth rate is measured as follows,
```math
OGT = 937 * F_{IVYWREL} − 335
```


## Setup

### Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

Getting the repository including sub-modules
```bash
git clone git@github.com:michoug/FertilityIslandMAGs.git
```

Create the main `snakemake` environment

```bash
# create venv
conda env create -f requirements.yaml -n "snakemake"
```

### Run Setup
* Place your preprocessed/trim reads (e.g. `sample_1_val_1.fq.gz` and `sample_2_val_2.fq.gz` files) in a `reads` folder
* Place the individual assemblies (e.g. `sample.fa`) into an `assembly` folder
* Modify the `config/config.yaml` file to change the different paths and eventually the different options
* Modify the `config/all_samples.txt` file to include your samples

### Without Slurm

`snakemake -s workflow/Snakefile --conda-prefix snakemake_envs --configfile config/config.yaml --cores 28 --use-conda -rp`

### With Slurm

This part was mainly taken from [@susheelbhanu](https://github.com/susheelbhanu/) [nomis_pipeline](https://github.com/susheelbhanu/nomis_pipeline)

* Modify the `slurm.yaml` file by checking `partition`, `qos` and `account` that heavily depends on your system
* Modify the `sbatch.sh` file by checking `#SBATCH -p`, `#SBATCH --qos=` and `#SBATCH -A` options that heavily depends on your system

`sbatch config/sbatch.sh`


##############################
# MODULES
import os, re
import glob
import pandas as pd
from snakemake.utils import validate

# ##############################
# # Parameters
# CORES=int(os.environ.get("CORES", 4))

##############################
# Config validation
validate(config, "../schemas/config.schema.yaml")

##############################
# default executable for snakemake
shell.executable("bash")


##############################
# working directory
workdir: config["work_dir"]


##############################
# Relevant directories
# ASSEMBLY_DIR=config["assembly_dir"]

WORK_DIR = config["work_dir"]
ASSEMBLY_DIR = config["assembly_dir"]
READS_DIR = config["reads_dir"]
RESULTS_DIR = config["results_dir"]

##############################
# Steps
STEPS = config["steps"]


##############################
# Input
with open("config/sample_list.txt") as _sample_list_fh:
    SAMPLES = [line.strip() for line in _sample_list_fh if line.strip()]

if not SAMPLES:
    raise ValueError("config/sample_list.txt contains no samples")

_duplicate_samples = {s for s in SAMPLES if SAMPLES.count(s) > 1}
if _duplicate_samples:
    raise ValueError(
        f"config/sample_list.txt contains duplicate sample names: {sorted(_duplicate_samples)}"
    )


##############################
# Memory per thread (MB). This is cluster-specific and can be overridden
# without editing config.yaml, e.g.:
#   snakemake ... --config memory_per_thread=8000
DEFAULT_MEMORY_PER_THREAD_MB = 7000


def get_mem_mb(wildcards, threads):
    return threads * config.get("memory_per_thread", DEFAULT_MEMORY_PER_THREAD_MB)

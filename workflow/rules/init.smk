##############################
# MODULES
import os, re
import glob
import pandas as pd

# ##############################
# # Parameters
# CORES=int(os.environ.get("CORES", 4))

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
SAMPLES = [line.strip() for line in open("config/sample_list.txt").readlines()]
# SAMPLES_2=SAMPLES


def get_mem_mb(wildcards, threads):
    return threads * config["memory_per_thread"]


##############################
# Global wildcard constraints
# {sample} must never contain "/": Snakemake's default wildcard regex (.+)
# matches it, which makes rule outputs ambiguous wherever a literal "/"
# follows {sample} elsewhere in the workflow - e.g.
# Bins/intermediate/{sample}/checkm2/... (checkm2.smk) vs.
# Bins/intermediate/{sample}/BinnerEval/{binner}/checkm2/... (binner_selection.smk),
# where the former's {sample} can absorb "<real_sample>/BinnerEval/<binner>"
# as one token. No sample/assembly identifier in this workflow legitimately
# contains "/".
#
# A narrower, rule-scoped `sample="[^/.]+"` constraint also lives on
# rule mappingHuman (trim_reads.smk) to resolve a similar ambiguity with
# removeHuman's {sample}.noHuman.bam output - scoped there rather than here
# because alreadyAssemble.smk derives ASSEMBLY (which feeds {sample}
# elsewhere) from arbitrary {assembly}.fna basenames, which may legitimately
# contain dots.
wildcard_constraints:
    sample="[^/]+",

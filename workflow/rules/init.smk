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
# GUNC database (downloaded in install_software_DB.smk, used by metawrap2.smk).
# Zenodo records and md5 sums are those of GUNC 1.1.1 (gunc/gunc_database.py).
GUNC_DOWNLOADS = {
    "progenomes_2.1": (
        "https://zenodo.org/records/19631736/files/gunc_db_progenomes2.1.dmnd.gz",
        "bc93a855e0760aad5c4e5f2d0e26da46",
        "447c9330056b02f29f30fe81fe4af4eb",
    ),
    "progenomes_3": (
        "https://zenodo.org/records/19631883/files/gunc_db_progenomes3.dmnd.gz",
        "a26b3be496017f291eece27740d0f37f",
        "b1a9347d219a632e9ecf4652dee9bdd1",
    ),
    "gtdb_95": (
        "https://zenodo.org/records/19631804/files/gunc_db_gtdb95.dmnd.gz",
        "14ed95a2fb1360925e28b7b55df14574",
        "3b502dc047efaafb9831a5eaec7617fd",
    ),
    "gtdb_214": (
        "https://zenodo.org/records/19632326/files/gunc_db_gtdb214.dmnd.gz",
        "3933007d83a7a85e4672295ee1f9d91f",
        "ac4d6304cab3d62a396703eeb039d3b7",
    ),
}
GUNC_CFG = config.get("gunc", {})
GUNC_ENABLED = GUNC_CFG.get("enabled", False)
GUNC_DB_NAME = GUNC_CFG.get("database", "progenomes_3")
if GUNC_ENABLED and GUNC_DB_NAME not in GUNC_DOWNLOADS:
    raise ValueError(
        f"gunc.database must be one of {sorted(GUNC_DOWNLOADS)}, got {GUNC_DB_NAME!r}"
    )
GUNC_URL, GUNC_GZ_MD5, GUNC_DMND_MD5 = GUNC_DOWNLOADS.get(
    GUNC_DB_NAME, GUNC_DOWNLOADS["progenomes_3"]
)
GUNC_DB = os.path.join(
    GUNC_CFG.get("db_dir", os.path.join(RESULTS_DIR, "DB/GUNC")),
    os.path.basename(GUNC_URL)[: -len(".gz")],
)


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

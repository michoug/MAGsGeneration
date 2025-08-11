# Pipeline for analysing the NOMIS MAGs
# Example call: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores 1 -rpn

##############################
# MODULES
import os, re
import glob
import shutil
import pandas as pd

##############################
# CONFIG
# can be overwritten by using --configfile <path to config> when calling snakemake
# configfile:"config/config.yaml"


include: "rules/init.smk"


ASSEMBLY = glob_wildcards(os.path.join(ASSEMBLY_DIR, "{assembly}.fna")).assembly

##############################
# TARGETS & RULES
# List of (main) targets to be created
TARGETS = []

if "install_software" in STEPS:

    include: "rules/install_software_DB.smk"

    TARGETS += ["status/install_software.done"]


if "trim_reads" in STEPS:

    include: "rules/trim_reads.smk"

    TARGETS += ["status/trimming.done"]


if "contigs_taxonomy" in STEPS:

    include: "rules/contigs_taxonomy_assembly.smk"

    TARGETS += ["status/taxonomy.done"]


if "check_assembly" in STEPS:

    include: "rules/filter_assembly_done.smk"

    TARGETS += ["status/checkAssembly.done"]


if "mags_generation" in STEPS:

    include: "rules/mags_generation.smk"

    TARGETS += ["status/mags_generation.done"]


if "binning_refinement" in STEPS:

    include: "rules/dastool.smk"

    TARGETS += ["status/binning_refinement.done"]


# CheckM
if "checkm" in STEPS:

    include: "rules/checkm2.smk"

    TARGETS += ["status/checkm2.done"]


# Refine

if "refine_all" in STEPS:

    include: "rules/refine.smk"

    TARGETS += ["status/refine.done"]


if "filter_checkm" in STEPS:

    include: "rules/filter_checkm2.smk"

    TARGETS += ["status/filter_checkm2.done"]


# Dereplication
if "dereplicate" in STEPS:

    include: "rules/galah.smk"

    TARGETS += ["status/dereplicate.done"]


if "annotations" in STEPS:

    include: "rules/annotations.smk"

    TARGETS += ["status/annotations.done"]


# Coverage
if "coverage" in STEPS:

    include: "rules/coverage.smk"

    TARGETS += ["status/coverage.done"]


# GTDBtk & CheckM2
if "quality" in STEPS:

    include: "rules/quality.smk"

    TARGETS += ["status/quality.done"]


# MG Thermometer
if "mgthermometer" in STEPS:

    include: "rules/mgthermometer.smk"

    TARGETS += ["status/mgthermometer.done"]


# EggNog-Mapper
if "eggnog" in STEPS:

    include: "rules/eggnog.smk"

    TARGETS += ["status/eggnog.done"]


if "singlem_all" in STEPS:

    include: "rules/singlem.smk"

    TARGETS += ["status/singlem.done"]


# No targets
if len(TARGETS) == 0:
    raise Exception("You are not serious. Nothing to be done? Really?")


rule NOMIS_MAGS:
    input:
        TARGETS,

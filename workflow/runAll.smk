# Pipeline for analysing MAGs
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


ASSEMBLY = SAMPLES

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


if "assembly_all" in STEPS:

    include: "rules/assembly.smk"

    TARGETS += ["status/assembly.done"]


if "contigs_taxonomy" in STEPS:

    include: "rules/contigs_taxonomy_mmseqs.smk"

    TARGETS += ["status/taxonomy.done"]


if "check_assembly" in STEPS:

    include: "rules/filter_assembly.smk"

    TARGETS += ["status/checkAssembly.done"]


# euk_mags_generation.smk reuses the prok rules ("use rule ... with:") and
# eukaryotes_quality.smk reuses the euk rules, so include the dependencies too
if (
    "prok_mags_generation" in STEPS
    or "euk_mags_generation" in STEPS
    or "eukaryotes_quality" in STEPS
):

    include: "rules/prok_mags_generation.smk"


if "euk_mags_generation" in STEPS or "eukaryotes_quality" in STEPS:

    include: "rules/euk_mags_generation.smk"


if "prok_mags_generation" in STEPS:
    TARGETS += ["status/prok_mags_generation.done"]


if "euk_mags_generation" in STEPS:
    TARGETS += ["status/euk_mags_generation.done"]


if "binner_selection" in STEPS:

    include: "rules/binner_selection.smk"

    TARGETS += ["status/binner_selection.done"]


if "binning_refinement" in STEPS:

    include: "rules/dastool.smk"

    TARGETS += ["status/binning_refinement.done"]


if "eukaryotes_quality" in STEPS:

    include: "rules/eukaryotes_quality.smk"

    TARGETS += ["status/eukaryotes_quality.done"]


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


if "syntracker" in STEPS:

    include: "rules/syntracker.smk"

    TARGETS += ["status/syntracker.done"]


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


if "report" in STEPS:

    include: "rules/report.smk"

    TARGETS += [os.path.join(RESULTS_DIR, "report/MAGs_summary.html")]


# No targets
if len(TARGETS) == 0:
    raise Exception("You are not serious. Nothing to be done? Really?")


rule MAGS_Generation:
    input:
        TARGETS,

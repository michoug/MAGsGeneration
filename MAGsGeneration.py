#!/usr/bin/env python3

import os
import shutil
import subprocess
import yaml
import argparse
import sys

# Optional steps
OPTIONAL_STEPS = {"syntracker", "singlem_all", "eggnog", "mgthermometer", "eukaryotes_filter"}

# Pipeline definitions
PIPELINES = {
    "runAll": "workflow/runAll.smk",
    "alreadyAssemble": "workflow/alreadyAssemble.smk"
}

# ----------------------------------------------------
# Install Pixi (if missing) + install all Pixi envs
# ----------------------------------------------------
def ensure_and_install_pixi_envs(pixi_project_path):
    pixi_toml = os.path.join(pixi_project_path, "pixi.toml")

    if not os.path.exists(pixi_toml):
        print(f"ℹ️ No pixi.toml found in {pixi_project_path}, skipping Pixi environment setup.")
        return

    # Ensure Pixi exists
    if shutil.which("pixi") is None:
        print("⚠️ Pixi not found. Attempting automatic installation…")
        try:
            subprocess.run(
                ["curl", "-fsSL", "https://pixi.sh/install.sh", "-o", "install_pixi.sh"],
                check=True
            )
            subprocess.run(["sh", "install_pixi.sh"], check=True)
            print("✅ Pixi installed successfully.\n")
        except Exception:
            print("❌ Failed to auto-install Pixi.\n"
                  "Please install manually: https://prefix.dev/pixi")
            sys.exit(1)

    # Install all Pixi environments
    print("📦 Installing Pixi environments defined in pixi.toml...")
    try:
        subprocess.run(["pixi", "install", "-a"], cwd=pixi_project_path, check=True)
        print("✅ Pixi environments installed.\n")
    except subprocess.CalledProcessError:
        print("❌ Pixi environment installation failed.")
        sys.exit(1)


# ---------------------------------------
# Helper functions
# ---------------------------------------
def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def save_config(config, path):
    with open(path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)


def run_snakemake(snakefile, configfile, cores=1, use_conda=False,
                  conda_prefix=None, dryrun=False, continue_run=False,
                  profile=None, extra_args=None):

    cmd = [
        "snakemake",
        "-s", snakefile,
        "--configfile", configfile,
        "--cores", str(cores)
    ]

    if use_conda:
        cmd += ["--use-conda"]
        if conda_prefix:
            cmd += ["--conda-prefix", conda_prefix]

    if dryrun:
        cmd.append("-n")

    if continue_run:
        cmd.append("--rerun-incomplete")

    if profile:
        cmd += ["--profile", profile]

    if extra_args:
        cmd += extra_args

    print(">>> Running Snakemake command:")
    print(" ".join(cmd), "\n")

    subprocess.run(cmd, check=True)


# -----------------------------------------
# Main
# -----------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Unified Snakemake Pipeline Launcher")

    parser.add_argument("-p", "--pipeline", choices=PIPELINES.keys(),
                        help="Pipeline to use (runAll, alreadyAssemble)")
    parser.add_argument("-c", "--config", required=True,
                        help="Main config YAML file")
    parser.add_argument("-o", "--output",
                        help="Output/results directory")
    parser.add_argument("--steps", nargs="+",
                        help="Optional steps to include")
    parser.add_argument("--use-conda", action="store_true",
                        help="Enable Snakemake conda environments")
    parser.add_argument("--conda-prefix",
                        default=os.path.join(os.environ.get("CONDA_PREFIX", "/tmp/conda"), "pipeline"),
                        help="Conda prefix for Snakemake")
    parser.add_argument("--cores", type=int, default=1,
                        help="Number of cores")
    parser.add_argument("--slurm", action="store_true",
                        help="Use Slurm profile")
    parser.add_argument("--continue", dest="continue_run", action="store_true",
                        help="Continue an existing run")
    parser.add_argument("-n", "--dryrun", action="store_true",
                        help="Dry run (Snakemake -n)")
    parser.add_argument("--snakemake-args", nargs=argparse.REMAINDER,
                        help="Extra arguments to pass directly to Snakemake (e.g., --unlock)")

    args = parser.parse_args()

    # Load config
    config = load_config(args.config)

    # ---------------------------------------------------
    # Determine pipeline
    # ---------------------------------------------------
    if args.continue_run:
        loaded_cfg = load_config(args.config)
        pipeline = loaded_cfg.get("pipeline")
        if not pipeline:
            raise ValueError("Existing config is missing 'pipeline'.")
    else:
        if not args.pipeline:
            raise ValueError("Pipeline must be specified when starting a new run.")
        pipeline = args.pipeline

    snakefile = PIPELINES[pipeline]

    # ---------------------------------------------------
    # Determine results/output directory
    # ---------------------------------------------------
    results_dir = (config.get("results_dir")
                   if args.continue_run
                   else os.path.abspath(args.output or config.get("results_dir")))

    if not results_dir:
        raise ValueError("No results directory provided (use -o PATH).")

    os.makedirs(results_dir, exist_ok=True)
    config["results_dir"] = results_dir

    # ---------------------------------------------------
    # Determine steps
    # ---------------------------------------------------
    all_steps = config.get("steps", [])
    mandatory = [s for s in all_steps if s not in OPTIONAL_STEPS]

    if args.steps:
        final_steps = mandatory + [s for s in args.steps if s in OPTIONAL_STEPS]
    else:
        final_steps = all_steps

    config["steps"] = final_steps
    config["pipeline"] = pipeline

    # Save run config
    run_config_path = os.path.join(results_dir, "config_run.yaml")
    save_config(config, run_config_path)

    # ---------------------------------------------------
    # Print step summary
    # ---------------------------------------------------
    optional_run = [s for s in final_steps if s in OPTIONAL_STEPS]
    mandatory_run = [s for s in final_steps if s not in OPTIONAL_STEPS]

    print("📌 Steps to run:")
    print("  Mandatory:", mandatory_run)
    print("  Optional :", optional_run, "\n")

    # -------------------------------------------------------------
    # Install Pixi if needed + install all pixi environments
    # -------------------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ensure_and_install_pixi_envs(script_dir)

    # Determine profile
    profile_path = os.path.join("profile") if args.slurm else None

    # ---------------------------------------------------
    # Run Snakemake
    # ---------------------------------------------------
    run_snakemake(
        snakefile=snakefile,
        configfile=run_config_path,
        cores=args.cores,
        use_conda=args.use_conda,
        conda_prefix=args.conda_prefix,
        dryrun=args.dryrun,
        continue_run=args.continue_run,
        profile=profile_path,
        extra_args=args.snakemake_args
    )


if __name__ == "__main__":
    main()

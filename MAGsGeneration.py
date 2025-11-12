import os
import shutil
import subprocess
import yaml
import argparse

# Optional steps
OPTIONAL_STEPS = {"syntracker", "singlem_all", "eggnog", "mgthermometer", "eukaryotes_filter"}
PIPELINES = {
    "runAll": "workflow/runAll.smk",
    "alreadyAssemble": "workflow/alreadyAssemble.smk"
}

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def save_config(config, path):
    with open(path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

def run_snakemake(snakefile, configfile, cores=1, use_conda=False, conda_prefix=None, dryrun=False, continue_run=False, profile=None):
    cmd = ["snakemake", "-s", snakefile, "--configfile", configfile, "--cores", str(cores)]
    if use_conda:
        cmd += ["--use-conda"]
        if conda_prefix:
            cmd += ["--conda-prefix", conda_prefix]
    if dryrun:
        cmd.append("--dry-run")
    if continue_run:
        cmd.append("--rerun-incomplete")
    if profile:
        cmd += ["--profile", profile]
    print(">>> Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pipeline", choices=PIPELINES.keys())
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--steps", nargs="+", help="Optional steps to include")
    parser.add_argument("-o", "--output", help="Results folder path")
    parser.add_argument("-n", "--dryrun", action="store_true", help="Perform a dry-run")
    parser.add_argument("--use-conda", action="store_true")
    parser.add_argument("--conda-prefix", default=os.path.join(os.environ.get("CONDA_PREFIX", "/tmp/conda"), "pipeline"))
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--slurm", action="store_true", help="Use Slurm profile for execution")
    parser.add_argument("--continue", dest="continue_run", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)

    # Determine pipeline
    pipeline = config.get("pipeline") if args.continue_run else args.pipeline
    if not pipeline:
        raise ValueError("Pipeline must be specified for new runs.")

    snakefile = PIPELINES[pipeline]

    # Determine results folder
    results_dir = config.get("results_dir") if args.continue_run else os.path.abspath(args.output or config.get("results_dir"))
    os.makedirs(results_dir, exist_ok=True)
    config["results_dir"] = results_dir

    # Filter steps
    all_steps = config.get("steps", [])
    mandatory = [s for s in all_steps if s not in OPTIONAL_STEPS]
    final_steps = mandatory + [s for s in args.steps if s in OPTIONAL_STEPS] if args.steps else all_steps
    config["steps"] = final_steps

    # Save run config in results folder
    run_config_path = os.path.join(results_dir, "config_run.yaml")
    save_config(config, run_config_path)

    # Print steps summary
    optional_run = [s for s in final_steps if s in OPTIONAL_STEPS]
    mandatory_run = [s for s in final_steps if s not in OPTIONAL_STEPS]
    print(f"Steps to run:\n  Mandatory: {mandatory_run}\n  Optional: {optional_run}")

    # Determine profile path for Snakemake
    profile_path = os.path.join("profile", "slurm") if args.slurm else None

    # Run Snakemake
    run_snakemake(
        snakefile=snakefile,
        configfile=run_config_path,
        cores=args.cores,
        use_conda=args.use_conda,
        conda_prefix=args.conda_prefix,
        dryrun=args.dryrun,
        continue_run=args.continue_run,
        profile=profile_path
    )

if __name__ == "__main__":
    main()

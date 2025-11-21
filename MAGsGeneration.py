import os
import subprocess
import yaml
import argparse
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Optional steps
OPTIONAL_STEPS = {
    "syntracker",
    "singlem_all",
    "eggnog",
    "mgthermometer",
    "eukaryotes_filter",
}
PIPELINES = {
    "runAll": "workflow/runAll.smk",
    "alreadyAssemble": "workflow/alreadyAssemble.smk",
}


def load_config(path):
    """Load and validate YAML configuration file."""
    try:
        with open(path) as f:
            config = yaml.safe_load(f)
        if not isinstance(config, dict):
            raise ValueError("Config must be a dictionary")
        return config
    except FileNotFoundError:
        raise ValueError(f"Config file not found: {path}")
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML in config: {e}")


def save_config(config, path):
    """Save configuration to YAML file."""
    try:
        with open(path, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)
        logger.info(f"Configuration saved to: {path}")
    except Exception as e:
        raise ValueError(f"Failed to save config to {path}: {e}")


def unlock_snakemake(snakefile, configfile):
    """Unlock Snakemake working directory before running."""
    cmd = ["snakemake", "-s", snakefile, "--configfile", configfile, "--unlock"]
    logger.info("Unlocking Snakemake pipeline...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Pipeline unlocked successfully")
    except subprocess.CalledProcessError as e:
        logger.warning(f"Unlock failed (may be okay if not locked): {e.stderr}")


def run_snakemake(
    snakefile,
    configfile,
    cores=1,
    use_conda=False,
    conda_prefix=None,
    dryrun=False,
    continue_run=False,
    profile=None,
):
    """Execute Snakemake pipeline with specified parameters."""
    cmd = ["snakemake", "-s", snakefile, "--configfile", configfile]

    # When using a profile (e.g., Slurm), cores should be set to "all"
    # The actual job resources are controlled by the profile configuration
    if profile:
        cmd += ["--cores", "all"]
        cmd += ["--profile", profile]
    else:
        cmd += ["--cores", str(cores)]

    if use_conda:
        cmd += ["--use-conda"]
        if conda_prefix:
            cmd += ["--conda-prefix", conda_prefix]

    if dryrun:
        cmd.append("--dry-run")

    if continue_run:
        cmd.append("--rerun-incomplete")

    logger.info(f"Running Snakemake: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        logger.info("Pipeline completed successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"Pipeline failed with exit code {e.returncode}")
        raise


def validate_config(config, continue_run=False):
    """Validate configuration has required keys."""
    if continue_run:
        return  # Skip validation for continuation runs

    required_keys = ["results_dir", "steps"]
    missing = [k for k in required_keys if k not in config]
    if missing:
        raise ValueError(f"Config missing required keys: {missing}")

    if not isinstance(config.get("steps"), list):
        raise ValueError("Config key 'steps' must be a list")


def validate_optional_steps(steps):
    """Validate that provided optional steps are valid."""
    if not steps:
        return

    # Allow special keywords
    if steps == ["all"] or steps == ["none"]:
        return

    invalid = set(steps) - OPTIONAL_STEPS
    if invalid:
        raise ValueError(
            f"Invalid optional steps: {invalid}. Valid options: {OPTIONAL_STEPS}"
        )


def resolve_optional_steps(args_steps, all_steps):
    """Resolve optional steps based on command line arguments.
    
    Returns the final list of optional steps to include.
    - If args_steps is None (not provided): default to 'all'
    - If args_steps is ['all']: include all optional steps
    - If args_steps is ['none']: include no optional steps
    - Otherwise: include specified optional steps
    """
    optional_in_config = [s for s in all_steps if s in OPTIONAL_STEPS]
    
    # Default to 'all' if --steps not provided
    if args_steps is None:
        logger.info("No --steps specified, defaulting to 'all' optional steps")
        return optional_in_config
    
    # Handle special keywords
    if args_steps == ["all"]:
        logger.info("Using 'all' optional steps")
        return optional_in_config
    
    if args_steps == ["none"]:
        logger.info("Using 'none' optional steps")
        return []
    
    # Filter to valid optional steps
    return [s for s in args_steps if s in OPTIONAL_STEPS]


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper for running Snakemake bioinformatics pipelines"
    )
    parser.add_argument(
        "-p",
        "--pipeline",
        choices=PIPELINES.keys(),
        help="Pipeline to run (required for new runs)",
    )
    parser.add_argument(
        "-c", "--config", required=True, help="Path to YAML configuration file"
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        help=f"Optional steps to include. Use 'all' for all optional steps (default), "
             f"'none' for no optional steps, or specify individual steps: {', '.join(OPTIONAL_STEPS)}",
    )
    parser.add_argument("-o", "--output", help="Results folder path")
    parser.add_argument(
        "-n",
        "--dryrun",
        action="store_true",
        help="Perform a dry-run without executing",
    )
    parser.add_argument(
        "--use-conda", action="store_true", help="Use conda environments"
    )
    parser.add_argument(
        "--conda-prefix",
        default=None,
        help="Conda prefix directory (default: $CONDA_PREFIX/pipeline or ./conda)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of cores for local execution (ignored when using --slurm, default: 1)",
    )
    parser.add_argument(
        "--slurm",
        action="store_true",
        help="Use Slurm profile for cluster execution (sets cores to 'all', job resources controlled by profile)",
    )
    parser.add_argument(
        "--profile",
        help="Path to Snakemake profile directory (overrides --slurm, sets cores to 'all')",
    )
    parser.add_argument(
        "--continue",
        dest="continue_run",
        action="store_true",
        help="Continue an incomplete run",
    )

    args = parser.parse_args()

    # Load configuration
    logger.info(f"Loading configuration from: {args.config}")
    config = load_config(args.config)
    validate_config(config, args.continue_run)

    # Determine pipeline
    if args.continue_run:
        pipeline = config.get("pipeline")
        if not pipeline:
            raise ValueError("Cannot continue: no pipeline specified in config")
        logger.info(f"Continuing pipeline: {pipeline}")
    else:
        if not args.pipeline:
            raise ValueError(
                "Pipeline must be specified for new runs using -p/--pipeline"
            )
        pipeline = args.pipeline
        config["pipeline"] = pipeline
        logger.info(f"Starting new pipeline: {pipeline}")

    if pipeline not in PIPELINES:
        raise ValueError(
            f"Unknown pipeline: {pipeline}. Choose from {list(PIPELINES.keys())}"
        )

    snakefile = PIPELINES[pipeline]

    # Determine results folder
    if args.continue_run:
        results_dir = config.get("results_dir")
        if not results_dir:
            raise ValueError("Cannot continue: no results_dir in config")
        results_dir = os.path.abspath(results_dir)
        if not os.path.exists(results_dir):
            raise ValueError(
                f"Cannot continue: results directory does not exist: {results_dir}"
            )
    else:
        results_dir = os.path.abspath(args.output or config.get("results_dir"))
        os.makedirs(results_dir, exist_ok=True)
        logger.info(f"Results directory: {results_dir}")

    config["results_dir"] = results_dir

    # Validate and resolve steps
    validate_optional_steps(args.steps)

    all_steps = config.get("steps", [])
    mandatory = [s for s in all_steps if s not in OPTIONAL_STEPS]
    
    # Resolve which optional steps to include (defaults to 'all')
    optional_to_include = resolve_optional_steps(args.steps, all_steps)
    
    final_steps = mandatory + optional_to_include
    config["steps"] = final_steps

    # Save run config in results folder
    run_config_path = os.path.join(results_dir, "config_run.yaml")

    # Check if continuing and config exists
    if args.continue_run:
        if not os.path.exists(run_config_path):
            raise ValueError(
                f"Cannot continue: no config_run.yaml found in {results_dir}"
            )
        logger.info(f"Using existing run config: {run_config_path}")
    else:
        save_config(config, run_config_path)

    # Print steps summary
    optional_run = [s for s in final_steps if s in OPTIONAL_STEPS]
    mandatory_run = [s for s in final_steps if s not in OPTIONAL_STEPS]
    logger.info(
        f"Steps to run:\n  Mandatory: {mandatory_run}\n  Optional: {optional_run}"
    )

    # Determine conda prefix
    if args.conda_prefix:
        conda_prefix = args.conda_prefix
    elif args.use_conda:
        conda_prefix = os.path.join(
            os.environ.get("CONDA_PREFIX", "./conda"), "pipeline"
        )
    else:
        conda_prefix = None

    # Determine profile path for Snakemake
    # Priority: --profile > --slurm > None
    profile_path = None
    if args.profile:
        profile_path = args.profile
        if not os.path.exists(profile_path):
            raise ValueError(f"Profile directory not found: {profile_path}")
        logger.info(f"Using custom profile: {profile_path}")
        logger.info("Note: --cores is ignored; job resources defined in profile config")
    elif args.slurm:
        profile_path = os.path.join("profile")
        if not os.path.exists(profile_path):
            raise ValueError(f"Slurm profile not found at: {profile_path}")
        logger.info(f"Using Slurm profile: {profile_path}")
        logger.info("Note: --cores is ignored; job resources defined in profile config")
    else:
        logger.info(f"Running locally with {args.cores} cores")

    # Always unlock before running
    unlock_snakemake(snakefile=snakefile, configfile=run_config_path)

    # Run Snakemake
    run_snakemake(
        snakefile=snakefile,
        configfile=run_config_path,
        cores=args.cores,
        use_conda=args.use_conda,
        conda_prefix=conda_prefix,
        dryrun=args.dryrun,
        continue_run=args.continue_run,
        profile=profile_path,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Error: {e}")
        exit(1)
"""Integration test: dry-run the real Snakemake DAG against a minimal fixture.

Unlike the other tests in this suite, these do not mock ``subprocess`` — they
invoke the real ``snakemake`` executable so that a broken rule graph (missing
input, unresolved wildcard, lint violation) is caught before it reaches a
production run. Skipped automatically if ``snakemake`` is not on PATH (e.g.
outside the pixi default environment).
"""

import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

snakemake_available = shutil.which("snakemake") is not None
pytestmark = pytest.mark.skipif(
    not snakemake_available, reason="snakemake is not on PATH"
)


def _write_fixture_config(tmp_path):
    """Build a minimal on-disk fixture and matching config.yaml for a dry-run."""
    work_dir = tmp_path / "work_dir"
    reads_dir = tmp_path / "reads"
    assembly_dir = tmp_path / "assembly"
    results_dir = tmp_path / "results"
    (work_dir / "config").mkdir(parents=True)
    reads_dir.mkdir()
    assembly_dir.mkdir()

    sample = "sampleA"
    (work_dir / "config" / "sample_list.txt").write_text(f"{sample}\n")
    (reads_dir / f"{sample}_R1_001.fastq.gz").write_bytes(b"")
    (reads_dir / f"{sample}_R2_001.fastq.gz").write_bytes(b"")

    config = {
        # trim_reads.smk references rules.downloadHuman (defined in
        # install_software_DB.smk) unconditionally, so that step must be
        # included whenever trim_reads is.
        "steps": ["install_software", "trim_reads", "assembly_all"],
        "work_dir": str(work_dir),
        "assembly_dir": str(assembly_dir),
        "results_dir": str(results_dir),
        "reads_dir": str(reads_dir),
        "memory_per_thread": 1000,
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    return config_path


def test_lint_does_not_error(tmp_path):
    """``snakemake --lint`` should not crash.

    ``--lint`` exits non-zero whenever it has *any* advisory finding, and this
    workflow has ~30 accepted "no conda: directive" findings (reproducibility
    is handled by Pixi instead — see PIPELINE_REVIEW.md), so asserting a
    zero exit code would never be satisfiable. What actually indicates a
    broken rule graph is a ``WorkflowError`` (e.g. a rule referencing another
    rule that isn't included for the given ``steps``), not an advisory lint
    finding — assert on that instead.
    """
    config_path = _write_fixture_config(tmp_path)
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/runAll.smk",
            "--configfile",
            str(config_path),
            "--lint",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert "WorkflowError" not in result.stderr, result.stdout + result.stderr


def test_dry_run_resolves_dag(tmp_path):
    """``snakemake -n`` should resolve the DAG without missing inputs."""
    config_path = _write_fixture_config(tmp_path)
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/runAll.smk",
            "--configfile",
            str(config_path),
            "-n",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

"""Tests verifying the Snakemake workflow files exist and pass snakefmt formatting."""

import os
import subprocess
import pytest
import glob

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKFLOW_DIR = os.path.join(REPO_ROOT, "workflow")
RULES_DIR = os.path.join(WORKFLOW_DIR, "rules")


def get_smk_files():
    """Return all .smk files in the workflow directory tree."""
    pattern = os.path.join(WORKFLOW_DIR, "**", "*.smk")
    return glob.glob(pattern, recursive=True)


def get_top_level_smk_files():
    """Return top-level workflow entry-point .smk files."""
    return glob.glob(os.path.join(WORKFLOW_DIR, "*.smk"))


# ---------------------------------------------------------------------------
# File existence tests
# ---------------------------------------------------------------------------


class TestWorkflowFilesExist:
    def test_runall_workflow_exists(self):
        assert os.path.exists(os.path.join(WORKFLOW_DIR, "runAll.smk"))

    def test_already_assemble_workflow_exists(self):
        assert os.path.exists(os.path.join(WORKFLOW_DIR, "alreadyAssemble.smk"))

    def test_rules_directory_exists(self):
        assert os.path.isdir(RULES_DIR)

    def test_init_rule_exists(self):
        assert os.path.exists(os.path.join(RULES_DIR, "init.smk"))

    def test_expected_rules_exist(self):
        expected_rules = [
            "init.smk",
            "trim_reads.smk",
            "assembly.smk",
            "filter_assembly.smk",
            "mags_generation.smk",
            "dastool.smk",
            "checkm2.smk",
            "refine.smk",
            "galah.smk",
            "coverage.smk",
            "quality.smk",
            "annotations.smk",
            "eggnog.smk",
            "singlem.smk",
            "syntracker.smk",
            "mgthermometer.smk",
            "eukaryotes_filter.smk",
        ]
        for rule_file in expected_rules:
            full_path = os.path.join(RULES_DIR, rule_file)
            assert os.path.exists(
                full_path
            ), f"Expected rule file not found: {rule_file}"

    def test_at_least_one_smk_file_found(self):
        smk_files = get_smk_files()
        assert len(smk_files) > 0, "No .smk files found in workflow directory"


# ---------------------------------------------------------------------------
# snakefmt format tests
# ---------------------------------------------------------------------------


def has_snakefmt():
    """Return True if snakefmt is available on PATH."""
    result = subprocess.run(["snakefmt", "--version"], capture_output=True, text=True)
    return result.returncode == 0


@pytest.mark.skipif(not has_snakefmt(), reason="snakefmt not installed")
class TestSnakefmtFormatting:
    @pytest.mark.parametrize("smk_file", get_smk_files())
    def test_smk_file_passes_snakefmt(self, smk_file):
        """Each .smk file should pass snakefmt --check (no formatting errors)."""
        result = subprocess.run(
            ["snakefmt", "--check", smk_file],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"snakefmt formatting issues in {os.path.relpath(smk_file, REPO_ROOT)}:\n"
            f"{result.stderr or result.stdout}"
        )

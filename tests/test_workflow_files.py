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


# ---------------------------------------------------------------------------
# Workflow structure tests
# ---------------------------------------------------------------------------


class TestWorkflowStructure:
    def test_workflow_directory_contains_runall(self):
        """Test that workflow directory contains runAll.smk."""
        runall = os.path.join(WORKFLOW_DIR, "runAll.smk")
        assert os.path.isfile(runall), "runAll.smk is not a file"

    def test_workflow_directory_contains_already_assemble(self):
        """Test that workflow directory contains alreadyAssemble.smk."""
        already = os.path.join(WORKFLOW_DIR, "alreadyAssemble.smk")
        assert os.path.isfile(already), "alreadyAssemble.smk is not a file"

    def test_rules_directory_is_subdirectory(self):
        """Test that rules is a subdirectory of workflow."""
        assert os.path.isdir(RULES_DIR), "rules is not a directory"
        assert RULES_DIR.startswith(WORKFLOW_DIR), "rules is not a subdirectory of workflow"


class TestRuleFileContent:
    def test_init_rule_has_setup_code(self):
        """Test that init.smk contains setup/initialization code."""
        init_file = os.path.join(RULES_DIR, "init.smk")
        with open(init_file, "r") as f:
            content = f.read()
            # init.smk contains global setup, not a rule definition
            assert len(content) > 0 and ("import" in content or "shell.executable" in content or "workdir" in content), \
                "init.smk does not contain setup code"

    def test_rule_files_contain_rule_definitions(self):
        """Test that rule files contain rule definitions."""
        sample_rules = ["trim_reads.smk", "assembly.smk", "coverage.smk"]
        
        for rule_file in sample_rules:
            full_path = os.path.join(RULES_DIR, rule_file)
            if os.path.exists(full_path):
                with open(full_path, "r") as f:
                    content = f.read()
                    assert "rule " in content, f"{rule_file} does not contain 'rule' keyword"


class TestWorkflowFileValidation:
    def test_all_rule_files_end_with_newline(self):
        """Test that all rule files end with a newline."""
        for smk_file in get_smk_files():
            with open(smk_file, "rb") as f:
                content = f.read()
                if len(content) > 0:
                    assert content.endswith(b"\n"), (
                        f"{smk_file} does not end with a newline"
                    )

    def test_rule_files_not_empty(self):
        """Test that all rule files have content."""
        for smk_file in get_smk_files():
            with open(smk_file, "r") as f:
                content = f.read().strip()
                assert len(content) > 0, f"{smk_file} is empty"

    def test_no_windows_line_endings(self):
        """Test that files use Unix line endings."""
        for smk_file in get_smk_files():
            with open(smk_file, "rb") as f:
                content = f.read()
                assert b"\r\n" not in content, (
                    f"{smk_file} has Windows line endings (CRLF)"
                )


class TestWorkflowIntegration:
    def test_both_main_workflows_exist(self):
        """Test that both main workflow entry points exist."""
        workflows = [
            os.path.join(WORKFLOW_DIR, "runAll.smk"),
            os.path.join(WORKFLOW_DIR, "alreadyAssemble.smk"),
        ]
        
        for workflow in workflows:
            assert os.path.exists(workflow), f"Main workflow missing: {workflow}"
            assert os.path.isfile(workflow), f"{workflow} is not a file"

    def test_rule_count_reasonable(self):
        """Test that a reasonable number of rules are defined."""
        smk_files = get_smk_files()
        assert len(smk_files) > 5, "Too few .smk rule files found"

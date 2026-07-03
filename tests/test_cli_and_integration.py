"""Integration tests for configuration file handling and pipelines.

This module includes both configuration validation tests and CLI integration tests
with mocker-based subprocess testing.
"""

import os
import sys
import pytest
import tempfile
import yaml
from unittest.mock import Mock, patch, MagicMock

# Add the project root to path so we can import MAGsGeneration
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from MAGsGeneration import load_config, PIPELINES, main


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_config_file(tmp_path):
    """Create a temporary config file for testing."""
    config = {
        "results_dir": str(tmp_path / "results"),
        "steps": ["trim_reads", "assembly_all", "mags_generation"],
        "work_dir": str(tmp_path / "work"),
        "assembly_dir": str(tmp_path / "assembly"),
        "reads_dir": str(tmp_path / "reads"),
        "memory_per_thread": 7000,
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f)
    return str(config_path)


@pytest.fixture
def tmp_profile_dir(tmp_path):
    """Create a temporary profile directory."""
    profile_dir = tmp_path / "profile"
    profile_dir.mkdir()
    # Create a minimal cluster config
    cluster_config = {
        "__default__": {
            "time": "00:30:00",
            "mem": 4000,
        }
    }
    cluster_path = profile_dir / "cluster.yaml"
    with open(cluster_path, "w") as f:
        yaml.safe_dump(cluster_config, f)
    return str(profile_dir)


# ---------------------------------------------------------------------------
# Pipeline and configuration validation tests
# ---------------------------------------------------------------------------


class TestPipelineDefinitions:
    def test_pipeline_choices_valid(self):
        """Test that pipeline choices match PIPELINES constant."""
        expected_pipelines = set(PIPELINES.keys())
        assert len(expected_pipelines) > 0
        assert "runAll" in expected_pipelines
        assert "alreadyAssemble" in expected_pipelines

    def test_pipeline_snakefiles_exist(self):
        """Test that all pipeline snakefiles exist."""
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for name, snakefile in PIPELINES.items():
            full_path = os.path.join(repo_root, snakefile)
            assert os.path.exists(full_path), f"Snakefile for pipeline '{name}' not found: {full_path}"
            assert os.path.isfile(full_path), f"Pipeline snakefile is not a file: {full_path}"


class TestConfigFileHandling:
    def test_load_and_validate_config(self, tmp_config_file):
        """Test loading a valid config file."""
        config = load_config(tmp_config_file)
        assert config is not None
        assert config["results_dir"]
        assert config["steps"]
        assert len(config["steps"]) > 0

    def test_config_has_required_fields(self, tmp_config_file):
        """Test that loaded config has required fields."""
        config = load_config(tmp_config_file)
        required = ["results_dir", "steps"]
        for field in required:
            assert field in config, f"Missing required field: {field}"

    def test_config_steps_is_list(self, tmp_config_file):
        """Test that steps field is a list."""
        config = load_config(tmp_config_file)
        assert isinstance(config["steps"], list), "steps must be a list"

    def test_config_results_dir_path_valid(self, tmp_config_file):
        """Test that results_dir path is a string."""
        config = load_config(tmp_config_file)
        assert isinstance(config["results_dir"], str)
        assert len(config["results_dir"]) > 0


# ---------------------------------------------------------------------------
# CLI argument parsing and workflow execution tests
# ---------------------------------------------------------------------------


class TestCLIArgumentParsing:
    def test_main_with_config_argument(self, mocker, tmp_config_file):
        """Test that main() accepts config argument."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        # Mock sys.argv
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "runAll"]
        mocker.patch.object(sys, "argv", test_args)
        
        # Call main (may exit, so we catch SystemExit)
        try:
            main()
        except SystemExit:
            pass

    def test_main_with_pipeline_selection(self, mocker, tmp_config_file):
        """Test that main() accepts pipeline selection."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "alreadyAssemble"]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass

    def test_main_with_dryrun_flag(self, mocker, tmp_config_file):
        """Test that main() accepts dryrun flag."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "runAll", "-n"]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass

    def test_main_with_cores_argument(self, mocker, tmp_config_file):
        """Test that main() accepts cores argument."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "runAll", "-j", "8"]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass

    def test_main_with_optional_steps(self, mocker, tmp_config_file):
        """Test that main() handles optional steps."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = [
            "MAGsGeneration.py",
            "-c",
            tmp_config_file,
            "-p",
            "runAll",
            "--optional-steps",
            "syntracker",
            "singlem_all",
        ]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass

    def test_main_calls_ensure_pixi(self, mocker, tmp_config_file):
        """Test that main() calls ensure_and_install_pixi_envs."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "runAll"]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass
        
        # Pixi setup should be called
        mock_ensure_pixi.assert_called()

    def test_main_calls_run_snakemake(self, mocker, tmp_config_file):
        """Test that main() calls run_snakemake."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "runAll"]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass
        
        # Snakemake should be called
        mock_run_snakemake.assert_called()


class TestWorkflowIntegration:
    def test_full_workflow_dry_run(self, mocker, tmp_config_file):
        """Test full workflow execution in dry-run mode."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = ["MAGsGeneration.py", "-c", tmp_config_file, "-p", "runAll", "-n"]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass
        
        # Both functions should be called
        mock_ensure_pixi.assert_called()
        mock_run_snakemake.assert_called()

    def test_workflow_with_profile(self, mocker, tmp_config_file, tmp_profile_dir):
        """Test workflow with custom profile."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = [
            "MAGsGeneration.py",
            "-c",
            tmp_config_file,
            "-p",
            "runAll",
            "--profile",
            tmp_profile_dir,
        ]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass
        
        # Should have called run_snakemake with profile
        mock_run_snakemake.assert_called()

    def test_workflow_with_conda(self, mocker, tmp_config_file):
        """Test workflow with conda environment."""
        mock_run_snakemake = mocker.patch("MAGsGeneration.run_snakemake")
        mock_ensure_pixi = mocker.patch("MAGsGeneration.ensure_and_install_pixi_envs")
        
        test_args = [
            "MAGsGeneration.py",
            "-c",
            tmp_config_file,
            "-p",
            "runAll",
            "--use-conda",
        ]
        mocker.patch.object(sys, "argv", test_args)
        
        try:
            main()
        except SystemExit:
            pass
        
        # Should have called run_snakemake
        mock_run_snakemake.assert_called()

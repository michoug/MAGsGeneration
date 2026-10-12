"""Unit tests for MAGsGeneration.py helper functions."""

import os
import sys
import pytest
import tempfile
import yaml

# Add the project root to path so we can import MAGsGeneration
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from MAGsGeneration import (
    load_config,
    save_config,
    validate_config,
    validate_optional_steps,
    resolve_optional_steps,
    OPTIONAL_STEPS,
    PIPELINES,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_config():
    return {
        "results_dir": "/tmp/results",
        "steps": ["trim_reads", "assembly_all", "prok_mags_generation"],
        "work_dir": "/tmp/work",
        "assembly_dir": "/tmp/assembly",
        "reads_dir": "/tmp/reads",
        "memory_per_thread": 7000,
    }


@pytest.fixture
def config_with_optional_steps(minimal_config):
    config = minimal_config.copy()
    config["steps"] = [
        "trim_reads",
        "assembly_all",
        "syntracker",
        "singlem_all",
        "eggnog",
    ]
    return config


@pytest.fixture
def tmp_yaml(tmp_path, minimal_config):
    """Write a minimal config to a temp YAML file and return its path."""
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(minimal_config))
    return str(path)


# ---------------------------------------------------------------------------
# load_config tests
# ---------------------------------------------------------------------------


class TestLoadConfig:
    def test_load_valid_config(self, tmp_yaml, minimal_config):
        config = load_config(tmp_yaml)
        assert config["results_dir"] == minimal_config["results_dir"]
        assert config["steps"] == minimal_config["steps"]

    def test_load_config_file_not_found(self):
        with pytest.raises(ValueError, match="Config file not found"):
            load_config("/nonexistent/path/config.yaml")

    def test_load_config_invalid_yaml(self, tmp_path):
        bad_yaml = tmp_path / "bad.yaml"
        bad_yaml.write_text("key: [unclosed bracket\n")
        with pytest.raises(ValueError, match="Invalid YAML"):
            load_config(str(bad_yaml))

    def test_load_config_not_a_dict(self, tmp_path):
        list_yaml = tmp_path / "list.yaml"
        list_yaml.write_text("- item1\n- item2\n")
        with pytest.raises(ValueError, match="Config must be a dictionary"):
            load_config(str(list_yaml))

    def test_load_config_returns_dict(self, tmp_yaml):
        config = load_config(tmp_yaml)
        assert isinstance(config, dict)


# ---------------------------------------------------------------------------
# save_config tests
# ---------------------------------------------------------------------------


class TestSaveConfig:
    def test_save_and_reload(self, tmp_path, minimal_config):
        out_path = str(tmp_path / "saved.yaml")
        save_config(minimal_config, out_path)
        assert os.path.exists(out_path)
        reloaded = load_config(out_path)
        assert reloaded == minimal_config

    def test_save_to_unwritable_path(self, minimal_config):
        with pytest.raises(ValueError, match="Failed to save config"):
            save_config(minimal_config, "/nonexistent_dir/config.yaml")

    def test_save_preserves_all_keys(self, tmp_path, minimal_config):
        out_path = str(tmp_path / "out.yaml")
        save_config(minimal_config, out_path)
        with open(out_path) as f:
            data = yaml.safe_load(f)
        assert set(data.keys()) == set(minimal_config.keys())


# ---------------------------------------------------------------------------
# validate_config tests
# ---------------------------------------------------------------------------


class TestValidateConfig:
    def test_valid_config_passes(self, minimal_config):
        validate_config(minimal_config)  # Should not raise

    def test_missing_results_dir_raises(self, minimal_config):
        del minimal_config["results_dir"]
        with pytest.raises(ValueError, match="results_dir"):
            validate_config(minimal_config)

    def test_missing_steps_raises(self, minimal_config):
        del minimal_config["steps"]
        with pytest.raises(ValueError, match="steps"):
            validate_config(minimal_config)

    def test_steps_not_a_list_raises(self, minimal_config):
        minimal_config["steps"] = "trim_reads"
        with pytest.raises(ValueError, match="must be a list"):
            validate_config(minimal_config)

    def test_continue_run_skips_validation(self):
        # When continue_run=True, even an empty config should be accepted
        validate_config({}, continue_run=True)  # Should not raise

    def test_both_keys_missing_lists_all(self):
        with pytest.raises(ValueError, match="results_dir") as exc_info:
            validate_config({})
        assert "steps" in str(exc_info.value)


# ---------------------------------------------------------------------------
# validate_optional_steps tests
# ---------------------------------------------------------------------------


class TestValidateOptionalSteps:
    def test_none_steps_passes(self):
        validate_optional_steps(None)  # Should not raise

    def test_empty_list_passes(self):
        validate_optional_steps([])  # Should not raise

    def test_all_keyword_passes(self):
        validate_optional_steps(["all"])  # Should not raise

    def test_none_keyword_passes(self):
        validate_optional_steps(["none"])  # Should not raise

    def test_valid_optional_steps_pass(self):
        for step in OPTIONAL_STEPS:
            validate_optional_steps([step])  # Should not raise

    def test_all_valid_optional_steps_pass(self):
        validate_optional_steps(list(OPTIONAL_STEPS))  # Should not raise

    def test_invalid_step_raises(self):
        with pytest.raises(ValueError, match="Invalid optional steps"):
            validate_optional_steps(["not_a_real_step"])

    def test_mix_valid_and_invalid_raises(self):
        valid = list(OPTIONAL_STEPS)[:1]
        with pytest.raises(ValueError, match="Invalid optional steps"):
            validate_optional_steps(valid + ["bogus_step"])


# ---------------------------------------------------------------------------
# resolve_optional_steps tests
# ---------------------------------------------------------------------------


class TestResolveOptionalSteps:
    def _steps_with_optional(self):
        return ["trim_reads", "assembly_all", "syntracker", "singlem_all", "eggnog"]

    def test_none_returns_all_optional_in_config(self):
        all_steps = self._steps_with_optional()
        result = resolve_optional_steps(None, all_steps)
        # Should return all optional steps found in config
        expected = [s for s in all_steps if s in OPTIONAL_STEPS]
        assert set(result) == set(expected)

    def test_all_keyword_returns_all_optional_in_config(self):
        all_steps = self._steps_with_optional()
        result = resolve_optional_steps(["all"], all_steps)
        expected = [s for s in all_steps if s in OPTIONAL_STEPS]
        assert set(result) == set(expected)

    def test_none_keyword_returns_empty(self):
        all_steps = self._steps_with_optional()
        result = resolve_optional_steps(["none"], all_steps)
        assert result == []

    def test_specific_steps_returns_only_those(self):
        all_steps = self._steps_with_optional()
        result = resolve_optional_steps(["syntracker"], all_steps)
        assert result == ["syntracker"]

    def test_non_optional_steps_filtered_out(self):
        all_steps = self._steps_with_optional()
        # trim_reads is mandatory, should be filtered out when passed as optional step
        result = resolve_optional_steps(["trim_reads"], all_steps)
        assert "trim_reads" not in result

    def test_no_optional_steps_in_config(self):
        mandatory_only = ["trim_reads", "assembly_all", "prok_mags_generation"]
        result = resolve_optional_steps(None, mandatory_only)
        assert result == []

    def test_default_none_equals_all(self):
        all_steps = self._steps_with_optional()
        result_none = resolve_optional_steps(None, all_steps)
        result_all = resolve_optional_steps(["all"], all_steps)
        assert set(result_none) == set(result_all)


# ---------------------------------------------------------------------------
# OPTIONAL_STEPS and PIPELINES constants
# ---------------------------------------------------------------------------


class TestConstants:
    def test_optional_steps_is_set(self):
        assert isinstance(OPTIONAL_STEPS, set)
        assert len(OPTIONAL_STEPS) > 0

    def test_expected_optional_steps_present(self):
        expected = {
            "syntracker",
            "singlem_all",
            "eggnog",
            "mgthermometer",
            "euk_mags_generation",
            "eukaryotes_quality",
        }
        assert expected.issubset(OPTIONAL_STEPS)

    def test_pipelines_is_dict(self):
        assert isinstance(PIPELINES, dict)
        assert len(PIPELINES) > 0

    def test_expected_pipelines_present(self):
        assert "runAll" in PIPELINES
        assert "alreadyAssemble" in PIPELINES

    def test_pipeline_snakefiles_exist(self):
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for name, snakefile in PIPELINES.items():
            full_path = os.path.join(repo_root, snakefile)
            assert os.path.exists(
                full_path
            ), f"Snakefile for pipeline '{name}' not found: {full_path}"


# ---------------------------------------------------------------------------
# unlock_snakemake and run_snakemake tests simplified (requires subprocess to work)
# ---------------------------------------------------------------------------

# These functions call subprocess.run; we test them by mocking subprocess.run via pytest-mock's `mocker` fixture.
# Detailed integration behavior is covered elsewhere; these tests focus on CLI argument-to-command translation.

# ---------------------------------------------------------------------------
# Error handling and edge cases
# ---------------------------------------------------------------------------


class TestConfigValidationEdgeCases:
    def test_config_with_extra_keys_passes(self, minimal_config):
        """Config with extra keys should pass validation."""
        minimal_config["extra_key"] = "extra_value"
        validate_config(minimal_config)  # Should not raise

    def test_empty_steps_list_passes(self, minimal_config):
        """Config with empty steps list should still be valid."""
        minimal_config["steps"] = []
        validate_config(minimal_config)  # Should not raise

    def test_steps_with_non_existent_steps_passes(self, minimal_config):
        """Config with non-existent steps should pass (validation is lenient)."""
        minimal_config["steps"] = ["trim_reads", "non_existent_step"]
        validate_config(minimal_config)  # Should not raise


class TestResolveOptionalStepsEdgeCases:
    def test_resolve_with_mixed_case_steps(self):
        """Test that step resolution works correctly with actual step names."""
        all_steps = ["trim_reads", "assembly_all", "syntracker", "eggnog"]
        result = resolve_optional_steps(["syntracker", "eggnog"], all_steps)
        assert set(result) == {"syntracker", "eggnog"}

    def test_resolve_filters_non_optional_correctly(self):
        """Test that mandatory steps are properly filtered."""
        all_steps = ["trim_reads", "assembly_all", "syntracker", "singlem_all"]
        result = resolve_optional_steps(["trim_reads", "syntracker"], all_steps)
        # trim_reads is mandatory, should not be in result
        assert "trim_reads" not in result
        assert "syntracker" in result

    def test_resolve_with_only_mandatory_steps(self):
        """Test resolution when only mandatory steps are provided."""
        all_steps = ["trim_reads", "assembly_all", "prok_mags_generation"]
        result = resolve_optional_steps(["trim_reads"], all_steps)
        assert result == []


# ---------------------------------------------------------------------------
# Subprocess-based functions with mocker
# ---------------------------------------------------------------------------


class TestUnlockSnakemake:
    def test_unlock_snakemake_calls_subprocess(self, mocker):
        """Test that unlock_snakemake calls subprocess.run with correct args."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import unlock_snakemake
        unlock_snakemake("test.smk", "test_config.yaml")
        
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "snakemake" in call_args
        assert "--unlock" in call_args

    def test_unlock_snakemake_with_specific_files(self, mocker):
        """Test that unlock_snakemake passes correct snakefile and config."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import unlock_snakemake
        snakefile = "workflow/runAll.smk"
        configfile = "config.yaml"
        unlock_snakemake(snakefile, configfile)
        
        call_args = mock_run.call_args[0][0]
        assert snakefile in call_args
        assert configfile in call_args


class TestRunSnakemake:
    def test_run_snakemake_calls_subprocess(self, mocker):
        """Test that run_snakemake calls subprocess.run."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            cores=4,
            dryrun=True,
            use_conda=True,
        )
        
        mock_run.assert_called_once()

    def test_run_snakemake_with_dryrun(self, mocker):
        """Test that dryrun flag is passed to snakemake."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            dryrun=True,
        )
        
        call_args = mock_run.call_args[0][0]
        assert "--dry-run" in call_args

    def test_run_snakemake_with_cores(self, mocker):
        """Test that cores parameter is passed."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            cores=8,
        )
        
        call_args = mock_run.call_args[0][0]
        assert "-j" in call_args or "--cores" in call_args or "8" in str(call_args)

    def test_run_snakemake_with_conda(self, mocker):
        """Test that use_conda flag enables conda."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            use_conda=True,
        )
        
        call_args = mock_run.call_args[0][0]
        assert "--use-conda" in call_args

    def test_run_snakemake_with_profile(self, mocker):
        """Test that profile parameter is passed."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            profile="slurm",
        )
        
        call_args = mock_run.call_args[0][0]
        assert "--profile" in call_args
        assert "slurm" in call_args

    def test_run_snakemake_with_continue_run(self, mocker):
        """Test that continue_run flag is handled."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            continue_run=True,
        )
        
        call_args = mock_run.call_args[0][0]
        assert "--rerun-incomplete" in call_args

    def test_run_snakemake_with_conda_prefix(self, mocker):
        """Test that conda_prefix is passed."""
        mock_run = mocker.patch("subprocess.run")
        
        from MAGsGeneration import run_snakemake
        run_snakemake(
            snakefile="test.smk",
            configfile="config.yaml",
            use_conda=True,
            conda_prefix="/tmp/conda",
        )
        
        call_args = mock_run.call_args[0][0]
        assert "--conda-prefix" in call_args
        assert "/tmp/conda" in call_args

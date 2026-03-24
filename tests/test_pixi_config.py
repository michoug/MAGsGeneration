"""Tests verifying the pixi.toml configuration structure."""

import os
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIXI_TOML = os.path.join(REPO_ROOT, "pixi.toml")


def load_pixi_toml():
    """Load and return the pixi.toml as a dict. Requires tomllib (Python 3.11+) or tomli."""
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            pytest.skip("tomllib/tomli not available – skipping TOML tests")
    with open(PIXI_TOML, "rb") as f:
        return tomllib.load(f)


# ---------------------------------------------------------------------------
# pixi.toml existence
# ---------------------------------------------------------------------------


class TestPixiTomlExists:
    def test_pixi_toml_file_exists(self):
        assert os.path.exists(PIXI_TOML), "pixi.toml not found at repository root"

    def test_pixi_lock_file_exists(self):
        lock_path = os.path.join(REPO_ROOT, "pixi.lock")
        assert os.path.exists(lock_path), "pixi.lock not found at repository root"


# ---------------------------------------------------------------------------
# pixi.toml structure
# ---------------------------------------------------------------------------


class TestPixiTomlStructure:
    @pytest.fixture(autouse=True)
    def pixi_config(self):
        self.config = load_pixi_toml()

    def test_workspace_section_exists(self):
        assert "workspace" in self.config

    def test_workspace_has_name(self):
        assert "name" in self.config["workspace"]

    def test_workspace_has_channels(self):
        assert "channels" in self.config["workspace"]
        assert len(self.config["workspace"]["channels"]) > 0

    def test_workspace_has_platforms(self):
        assert "platforms" in self.config["workspace"]
        assert "linux-64" in self.config["workspace"]["platforms"]

    def test_feature_section_exists(self):
        assert "feature" in self.config

    def test_environments_section_exists(self):
        assert "environments" in self.config

    def test_each_environment_has_feature(self):
        environments = self.config.get("environments", {})
        features = self.config.get("feature", {})
        for env_name, env_def in environments.items():
            feature_list = env_def.get("features", [])
            for feat in feature_list:
                assert feat in features, (
                    f"Environment '{env_name}' references feature '{feat}' "
                    f"which is not defined in [feature]"
                )

    def test_expected_environments_exist(self):
        environments = self.config.get("environments", {})
        expected_envs = [
            "megahit",
            "metabat",
            "coverm",
            "concoct",
            "semibin",
            "vamb",
            "dastool",
            "checkm2",
            "galah",
            "gtdbtk",
            "eggnog",
            "annotations",
            "singlem",
            "fastp",
            "trimgalore",
            "busco",
            "rosella",
        ]
        for env in expected_envs:
            assert (
                env in environments
            ), f"Expected environment '{env}' not found in pixi.toml"

    def test_snakemake_in_dependencies(self):
        deps = self.config.get("dependencies", {})
        assert "snakemake" in deps, "snakemake should be in top-level [dependencies]"

    def test_snakefmt_in_dependencies(self):
        deps = self.config.get("dependencies", {})
        assert "snakefmt" in deps, "snakefmt should be in top-level [dependencies]"

    def test_each_feature_has_dependencies(self):
        features = self.config.get("feature", {})
        for feat_name, feat_def in features.items():
            assert (
                "dependencies" in feat_def
            ), f"Feature '{feat_name}' has no 'dependencies' key"
            assert (
                len(feat_def["dependencies"]) > 0
            ), f"Feature '{feat_name}' has an empty 'dependencies' dict"

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


# ---------------------------------------------------------------------------
# Additional structure validation tests
# ---------------------------------------------------------------------------


class TestPixiTomlDependencies:
    @pytest.fixture(autouse=True)
    def pixi_config(self):
        self.config = load_pixi_toml()

    def test_top_level_dependencies_section_exists(self):
        """Test that top-level dependencies section exists."""
        assert "dependencies" in self.config, "No top-level [dependencies] section found"
        assert len(self.config["dependencies"]) > 0, "Top-level dependencies section is empty"

    def test_essential_tools_in_dependencies(self):
        """Test that essential tools are in dependencies."""
        deps = self.config.get("dependencies", {})
        essential = ["snakemake"]
        for tool in essential:
            assert tool in deps, f"Essential tool '{tool}' not in top-level dependencies"


class TestPixiFeatures:
    @pytest.fixture(autouse=True)
    def pixi_config(self):
        self.config = load_pixi_toml()

    def test_feature_names_valid(self):
        """Test that feature names are valid identifiers."""
        features = self.config.get("feature", {})
        for feat_name in features.keys():
            assert feat_name.isidentifier() or "_" in feat_name, (
                f"Feature name '{feat_name}' is not a valid identifier"
            )

    def test_feature_dependencies_are_dicts(self):
        """Test that feature dependencies are dictionaries."""
        features = self.config.get("feature", {})
        for feat_name, feat_def in features.items():
            deps = feat_def.get("dependencies")
            assert isinstance(deps, dict), (
                f"Feature '{feat_name}' dependencies are not a dict"
            )

    def test_no_duplicate_features(self):
        """Test that there are no duplicate feature names (case-insensitive)."""
        features = self.config.get("feature", {})
        feature_lower = [name.lower() for name in features.keys()]
        assert len(feature_lower) == len(set(feature_lower)), (
            "Duplicate feature names found (case-insensitive)"
        )


class TestPixiEnvironments:
    @pytest.fixture(autouse=True)
    def pixi_config(self):
        self.config = load_pixi_toml()

    def test_environment_definitions_are_dicts(self):
        """Test that environment definitions are dictionaries."""
        environments = self.config.get("environments", {})
        for env_name, env_def in environments.items():
            assert isinstance(env_def, dict), (
                f"Environment '{env_name}' definition is not a dict"
            )

    def test_environments_have_features_or_dependencies(self):
        """Test that each environment has either features or dependencies."""
        environments = self.config.get("environments", {})
        for env_name, env_def in environments.items():
            has_features = "features" in env_def and len(env_def["features"]) > 0
            has_dependencies = "dependencies" in env_def and len(env_def["dependencies"]) > 0
            assert has_features or has_dependencies, (
                f"Environment '{env_name}' has neither features nor dependencies"
            )

    def test_environment_features_exist(self):
        """Test that all referenced features actually exist."""
        environments = self.config.get("environments", {})
        features = self.config.get("feature", {})
        
        for env_name, env_def in environments.items():
            feature_list = env_def.get("features", [])
            for feat in feature_list:
                assert feat in features, (
                    f"Environment '{env_name}' references undefined feature '{feat}'"
                )

    def test_no_duplicate_environments(self):
        """Test that there are no duplicate environment names."""
        environments = self.config.get("environments", {})
        env_names = list(environments.keys())
        assert len(env_names) == len(set(env_names)), (
            "Duplicate environment names found"
        )


class TestPixiTomlValidation:
    @pytest.fixture(autouse=True)
    def pixi_config(self):
        self.config = load_pixi_toml()

    def test_workspace_name_not_empty(self):
        """Test that workspace name is not empty."""
        workspace = self.config.get("workspace", {})
        name = workspace.get("name", "")
        assert name, "Workspace name is empty"
        assert len(name) > 0, "Workspace name should not be empty"

    def test_channels_are_non_empty(self):
        """Test that channels list is non-empty."""
        workspace = self.config.get("workspace", {})
        channels = workspace.get("channels", [])
        assert len(channels) > 0, "Channels list is empty"

    def test_valid_channel_names(self):
        """Test that channel names look valid."""
        workspace = self.config.get("workspace", {})
        channels = workspace.get("channels", [])
        for channel in channels:
            # Channels should not be empty strings
            assert channel and len(channel) > 0, f"Invalid channel: {channel}"
            # Common channels
            valid_start = (
                channel.startswith("conda-forge") or
                channel.startswith("bioconda") or
                channel.startswith("http") or
                channel.startswith("file://") or
                "/" in channel
            )
            assert valid_start or channel in ["defaults"], (
                f"Channel name looks invalid: {channel}"
            )

    def test_platforms_contains_linux64(self):
        """Test that linux-64 platform is included."""
        workspace = self.config.get("workspace", {})
        platforms = workspace.get("platforms", [])
        assert "linux-64" in platforms, "linux-64 platform is missing"

    def test_no_conflicting_dependency_versions(self):
        """Test that dependency version constraints are valid."""
        deps = self.config.get("dependencies", {})
        for dep_name, dep_spec in deps.items():
            # dep_spec might be a string like "snakemake=7.32" or "python>=3.9"
            if isinstance(dep_spec, str):
                # Check for obviously invalid patterns
                assert not dep_spec.startswith("="), (
                    f"Dependency '{dep_name}' has invalid version spec: {dep_spec}"
                )
                # Should contain valid operator characters
                if any(op in dep_spec for op in ["=", ">", "<", "!"]):
                    pass  # Has version constraint - looks valid

    def test_pyproject_table_valid_if_exists(self):
        """Test that pyproject table is valid if it exists."""
        if "pyproject" in self.config:
            pyproject = self.config["pyproject"]
            # Should be a dict
            assert isinstance(pyproject, dict), "pyproject table is not a dict"


# ---------------------------------------------------------------------------
# pixi.lock validation tests
# ---------------------------------------------------------------------------


class TestPixiLockFile:
    def test_pixi_lock_is_readable(self):
        lock_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "pixi.lock"
        )
        try:
            with open(lock_path, "r") as f:
                content = f.read()
        except FileNotFoundError as e:
            pytest.skip(f"Could not read pixi.lock: {e}")
        assert len(content) > 0, "pixi.lock is empty"

    def test_pixi_lock_appears_valid(self):
        """Test that pixi.lock has expected structure."""
        lock_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "pixi.lock"
        )
        yaml = pytest.importorskip("yaml")
        try:
            with open(lock_path, "r", encoding="utf-8") as f:
                lock_data = yaml.safe_load(f)
        except FileNotFoundError as e:
            pytest.skip(f"Could not read pixi.lock: {e}")

        assert isinstance(lock_data, dict) and "version" in lock_data and "environments" in lock_data, (
            "pixi.lock appears invalid or missing expected keys"
        )

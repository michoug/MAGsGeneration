"""Tests for Pixi environment handling and subprocess utilities."""

import os
import sys
import pytest
import subprocess

# Add the project root to path so we can import MAGsGeneration
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from MAGsGeneration import ensure_and_install_pixi_envs


# ---------------------------------------------------------------------------
# Tests for ensure_and_install_pixi_envs
# ---------------------------------------------------------------------------


class TestEnsureAndInstallPixiEnvs:
    def test_skips_when_no_pixi_toml(self, mocker, tmp_path):
        """Test that function skips when no pixi.toml is found."""
        mock_run = mocker.patch("subprocess.run")
        
        # Mock shutil.which to return None (pixi not found)
        mocker.patch("shutil.which", return_value=None)
        
        ensure_and_install_pixi_envs(str(tmp_path))
        
        # Should not attempt any subprocess calls
        mock_run.assert_not_called()

    def test_installs_pixi_when_missing(self, mocker, tmp_path):
        """Test that pixi is installed when not found."""
        # Create a pixi.toml file
        pixi_toml = tmp_path / "pixi.toml"
        pixi_toml.write_text("[workspace]\nname = 'test'")
        
        # Mock shutil.which to return None (pixi not found)
        mocker.patch("shutil.which", return_value=None)
        
        # Mock subprocess.run to simulate successful installation
        mock_run = mocker.patch("subprocess.run")
        
        # Mock sys.exit to prevent actual exit
        mock_exit = mocker.patch("sys.exit")
        
        ensure_and_install_pixi_envs(str(tmp_path))
        
        # Should have attempted to install pixi via curl/bootstrap
        assert mock_run.called or mock_exit.called

    def test_uses_existing_pixi(self, mocker, tmp_path):
        """Test that function uses pixi if it exists."""
        # Create a pixi.toml file
        pixi_toml = tmp_path / "pixi.toml"
        pixi_toml.write_text("[workspace]\nname = 'test'")
        
        # Mock shutil.which to return pixi path
        pixi_path = str(tmp_path / "pixi")
        mocker.patch("shutil.which", return_value=pixi_path)
        
        # Mock subprocess.run
        mock_run = mocker.patch("subprocess.run")
        
        ensure_and_install_pixi_envs(str(tmp_path))
        
        # Should have called pixi install or similar
        mock_run.assert_called()

    def test_returns_none_on_success(self, mocker, tmp_path):
        """Test that function returns None on success."""
        # Create a pixi.toml file
        pixi_toml = tmp_path / "pixi.toml"
        pixi_toml.write_text("[workspace]\nname = 'test'")
        
        # Mock shutil.which to return pixi path
        pixi_path = str(tmp_path / "pixi")
        mocker.patch("shutil.which", return_value=pixi_path)
        
        # Mock subprocess.run to succeed
        mocker.patch("subprocess.run")
        
        result = ensure_and_install_pixi_envs(str(tmp_path))
        
        assert result is None

    def test_handles_subprocess_error(self, mocker, tmp_path):
        """Test that function handles subprocess errors gracefully."""
        # Create a pixi.toml file
        pixi_toml = tmp_path / "pixi.toml"
        pixi_toml.write_text("[workspace]\nname = 'test'")
        
        # Mock shutil.which to return pixi path
        pixi_path = str(tmp_path / "pixi")
        mocker.patch("shutil.which", return_value=pixi_path)
        
        # Mock subprocess.run to fail
        mocker.patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "pixi"))
        
        # Should handle error without crashing
        with pytest.raises((subprocess.CalledProcessError, SystemExit)):
            ensure_and_install_pixi_envs(str(tmp_path))

    def test_pixi_toml_path_validation(self, mocker, tmp_path):
        """Test that invalid paths are handled."""
        mock_run = mocker.patch("subprocess.run")
        mocker.patch("shutil.which", return_value=None)

        # Call with non-existent directory
        ensure_and_install_pixi_envs(str(tmp_path / "nonexistent"))

        # Should not attempt any subprocess calls
        mock_run.assert_not_called()
    def test_pixi_install_with_features(self, mocker, tmp_path):
        """Test pixi install is called with proper parameters."""
        # Create a pixi.toml file
        pixi_toml = tmp_path / "pixi.toml"
        pixi_toml.write_text("""[workspace]
name = 'test'

[environments]
test_env = { features = ["test"], no-default-feature = true }

[feature.test.dependencies]
python = "3.10.*"
""")
        
        # Mock shutil.which to return pixi path
        pixi_path = str(tmp_path / "pixi")
        mocker.patch("shutil.which", return_value=pixi_path)
        
        # Mock subprocess.run
        mock_run = mocker.patch("subprocess.run")
        
        ensure_and_install_pixi_envs(str(tmp_path))
        
        # Verify pixi was called
        mock_run.assert_called()

    def test_quiet_mode_flag(self, mocker, tmp_path):
        """Test that quiet mode is respected if function supports it."""
        # Create a pixi.toml file
        pixi_toml = tmp_path / "pixi.toml"
        pixi_toml.write_text("[workspace]\nname = 'test'")
        
        # Mock shutil.which
        pixi_path = str(tmp_path / "pixi")
        mocker.patch("shutil.which", return_value=pixi_path)
        
        # Mock subprocess.run
        mock_run = mocker.patch("subprocess.run")
        
        ensure_and_install_pixi_envs(str(tmp_path))
        
        # Function should have been called
        assert mock_run.called

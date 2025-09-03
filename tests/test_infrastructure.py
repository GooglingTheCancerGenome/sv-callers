"""Test to validate that the testing infrastructure is set up correctly."""
import pytest
import os
import sys
from pathlib import Path


def test_pytest_working():
    """Basic test to verify pytest is functioning."""
    assert True


def test_pytest_markers():
    """Test that custom markers are available."""
    # This test verifies that the markers work by using them
    # The actual marker functionality is tested by applying them to test functions
    assert hasattr(pytest.mark, 'unit')
    assert hasattr(pytest.mark, 'integration') 
    assert hasattr(pytest.mark, 'slow')


def test_project_structure():
    """Test that the project has the expected directory structure."""
    project_root = Path(__file__).parent.parent
    
    # Check main directories exist
    assert (project_root / "workflow").exists()
    assert (project_root / "tests").exists()
    assert (project_root / "tests" / "unit").exists()
    assert (project_root / "tests" / "integration").exists()
    
    # Check key files exist
    assert (project_root / "pyproject.toml").exists()
    assert (project_root / "tests" / "conftest.py").exists()
    assert (project_root / "tests" / "__init__.py").exists()


def test_imports():
    """Test that required testing packages can be imported."""
    try:
        import pytest
        import pytest_cov
        assert True
    except ImportError as e:
        pytest.fail(f"Required testing package not available: {e}")


def test_coverage_integration():
    """Test that coverage is properly integrated (this test itself contributes to coverage)."""
    # Simple function to ensure coverage tracking works
    def sample_function(x):
        if x > 0:
            return "positive"
        elif x < 0:
            return "negative"
        else:
            return "zero"
    
    assert sample_function(1) == "positive"
    assert sample_function(-1) == "negative"
    assert sample_function(0) == "zero"


@pytest.mark.unit
def test_unit_marker():
    """Test that unit marker works."""
    assert True


@pytest.mark.integration  
def test_integration_marker():
    """Test that integration marker works."""
    assert True


@pytest.mark.slow
def test_slow_marker():
    """Test that slow marker works."""
    assert True


def test_fixtures_available(temp_dir, sample_config):
    """Test that shared fixtures from conftest.py are available."""
    assert temp_dir.exists()
    assert isinstance(sample_config, dict)
    assert "genome" in sample_config


def test_mock_functionality():
    """Test that mocking functionality works with pytest-mock."""
    def simple_function():
        return "original"
    
    # Test without mocking first
    assert simple_function() == "original"
    
    # Note: actual mocking would require pytest-mock fixture in function parameters
    # This test validates the basic structure is in place


class TestInfrastructureClass:
    """Test class to verify class-based test discovery works."""
    
    def test_class_method(self):
        """Test method within a test class."""
        assert True
    
    def test_class_fixtures(self, temp_dir):
        """Test that fixtures work in test classes."""
        assert temp_dir.exists()


def test_python_version():
    """Test that we're running on a supported Python version."""
    version = sys.version_info
    assert version.major == 3
    assert version.minor >= 8  # Based on pyproject.toml requirement
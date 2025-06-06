import pytest
import os
import sys

# Add the project root directory to Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Add any common fixtures here
@pytest.fixture
def sample_data():
    """Sample fixture that can be used across test files."""
    return {
        "lat": 10.77532775390349,
        "lon": 106.70647829173976
    } 
import pytest
from vgrid.conversion.latlon2dggs import latlon2h3

def test_latlon2h3_basic(sample_data):
    """Test basic lat/lon to H3 conversion."""
    lat, lon = sample_data["lat"], sample_data["lon"]
    h3_index = latlon2h3(lat, lon, res=9)
    assert h3_index is not None
    assert isinstance(h3_index, str)
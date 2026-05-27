"""DGGS resampling utilities (area-weighted and nearest-neighbour transfer)."""

from vgrid.conversion.dggsresample.dggsresample import (
    generate_grid,
    get_nearest_resolution,
    dggsresample,
    dggsresample_cli,
    resampling,
)
from vgrid.utils.io import process_input_data_resample

__all__ = [
    "generate_grid",
    "get_nearest_resolution",
    "dggsresample",
    "dggsresample_cli",
    "resampling",
    "process_input_data_resample",
]

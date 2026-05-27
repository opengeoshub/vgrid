"""
This module provides functions for generating statistics for MGRS DGGS cells.
"""

import pandas as pd
import argparse
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import characteristic_length_scale

min_res = DGGS_TYPES["mgrs"]["min_res"]
max_res = DGGS_TYPES["mgrs"]["max_res"]


def mgrs_metrics(resolution, unit: str = "m"):
    """
    Calculate metrics for MGRS DGGS cells.

    Args:
        resolution: Resolution level (0-5)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, avg_edge_len_in_unit, avg_cell_area_in_unit_squared)
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    latitude_degrees = 8  # The latitude span of each GZD cell in degrees
    longitude_degrees = 6  # The longitude span of each GZD cell in degrees
    km_per_degree = 111  # Approximate kilometers per degree of latitude/longitude
    gzd_cells = 1200  # Total number of GZD cells

    # Convert degrees to kilometers
    latitude_span = latitude_degrees * km_per_degree
    longitude_span = longitude_degrees * km_per_degree

    # Calculate cell size in kilometers based on resolution
    # Resolution 1: 100 km, each subsequent resolution divides by 10
    cell_size_km = 100 / (10 ** (resolution))   
    # Calculate number of cells in latitude and longitude for the chosen cell size
    cells_latitude = latitude_span / cell_size_km
    cells_longitude = longitude_span / cell_size_km

    # Total number of cells for each GZD cell
    cells_per_gzd_cell = cells_latitude * cells_longitude

    # Total number of cells for all GZD cells
    num_cells = cells_per_gzd_cell * gzd_cells
    avg_edge_len = cell_size_km  # in km
    avg_cell_area = avg_edge_len**2  # in km2
    cls = characteristic_length_scale(
        avg_cell_area * (10 ** (-6)), unit=unit
    )  # convert avg_cell_area to m2 before calling characteristic_length_scale

    if unit == "m":
        avg_edge_len = cell_size_km * (10**3)  # Convert km to m
        avg_cell_area = avg_cell_area * (10**6)

    return num_cells, avg_edge_len, avg_cell_area, cls


def mgrsstats(unit: str = "m"):
    """
    Generate statistics for MGRS DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing MGRS DGGS statistics with columns:
            - resolution: Resolution level (0-5)
            - number_of_cells: Number of cells at each resolution
            - avg_edge_len_{unit}: Average edge length in the given unit
            - avg_cell_area_{unit}2: Average cell area in the squared unit
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # Initialize lists to store data
    resolutions = []
    num_cells_list = []
    avg_edge_lens = []
    avg_cell_areas = []
    cls_list = []
    for resolution in range(min_res, max_res + 1):
        num_cells, avg_edge_len, avg_cell_area, cls = mgrs_metrics(resolution, unit=unit)
        resolutions.append(resolution)
        num_cells_list.append(num_cells)
        avg_edge_lens.append(avg_edge_len)
        avg_cell_areas.append(avg_cell_area)
        cls_list.append(cls)
    # Create DataFrame
    # Build column labels with unit awareness (lower case)
    avg_edge_len = f"avg_edge_len_{unit}"
    unit_area_label = {"m": "m2", "km": "km2"}[unit]
    avg_cell_area = f"avg_cell_area_{unit_area_label}"
    cls_label = f"cls_{unit}"
    df = pd.DataFrame(
        {
            "resolution": resolutions,
            "number_of_cells": num_cells_list,
            avg_edge_len: avg_edge_lens,
            avg_cell_area: avg_cell_areas,
            cls_label: cls_list,
        }
    )

    return df


def mgrsstats_cli():
    """
    Command-line interface for generating MGRS DGGS statistics.

    CLI options:
      -unit, --unit {m,km}
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()  # type: ignore

    unit = args.unit

    print("Resolution 0: 100 x 100 km")
    print("Resolution 1: 10 x 10 km")
    print("2 <= Resolution <= 5 = Finer subdivisions (1 x 1 km, 0.1 x 0.11 km, etc.)")

    # Get the DataFrame
    df = mgrsstats(unit=unit)

    # Display the DataFrame
    print(df)


if __name__ == "__main__":
    mgrsstats_cli()

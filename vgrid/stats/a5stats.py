"""
T`h`is module provides functions for generating statistics for A5 DGGS cells.
"""

import math
import json
import pandas as pd
import numpy as np
import argparse
import geopandas as gpd
from a5.core.cell_info import get_num_cells, cell_area
from vgrid.generator.a5grid import a5grid
from vgrid.utils.geometry import (
    check_crossing_geom,
    characteristic_length_scale,
    convexhull_from_lambert,
    get_area_perimeter_from_lambert,
    get_cells_area,
)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm
from vgrid.utils.constants import DGGS_TYPES, VMIN_PEN, VMAX_PEN, VCENTER_PEN

min_res = DGGS_TYPES["a5"]["min_res"]
max_res = DGGS_TYPES["a5"]["max_res"]


def a5_metrics(resolution: int, unit: str = "m"):  # length unit is m, area unit is m2
    """
    Calculate metrics for A5 DGGS cells at a given resolution.

    Args:
        resolution: Resolution level (0-29)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, edge_length_in_unit, cell_area_in_unit_squared)
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    num_cells = get_num_cells(resolution)
    avg_cell_area = cell_area(resolution)  # cell_area returns area in m²

    # Calculate edge length in meters
    k = math.sqrt(5 * (5 + 2 * math.sqrt(5)))
    avg_edge_len = math.sqrt(4 * avg_cell_area / k)  # edge length in m
    cls = characteristic_length_scale(avg_cell_area, unit=unit)
    # Convert to requested unit
    if unit == "km":
        avg_edge_len = avg_edge_len / (10**3)
        avg_cell_area = avg_cell_area / (10**6)

    return num_cells, avg_edge_len, avg_cell_area, cls


def a5stats(unit: str = "m"):
    """
    Generate statistics for A5 DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing A5 DGGS statistics with columns:
            - Resolution: Resolution level (0-29)
            - Number_of_Cells: Number of cells at each resolution
            - Avg_Edge_Length_{unit}: Average edge length in the given unit
            - CLS: Characteristic length scale in the given unit
            - Avg_Cell_Area_{unit}2: Average cell area in the squared unit
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # Derive bounds from central constants registry

    # Initialize lists to store data
    resolutions = []
    num_cells_list = []
    avg_edge_lens = []
    avg_cell_areas = []
    cls_list = []
    for res in range(min_res, max_res + 1):
        num_cells, avg_edge_len, avg_cell_area, cls = a5_metrics(
            res, unit=unit
        )  # length unit is m, area unit is m2
        resolutions.append(res)
        num_cells_list.append(num_cells)
        avg_edge_lens.append(avg_edge_len)
        avg_cell_areas.append(avg_cell_area)
        cls_list.append(cls)

    # Create DataFrame
    # Build column labels with unit awareness
    avg_edge_len = f"avg_edge_len_{unit}"
    unit_area_label = {"m": "m2", "km": "km2"}[unit]
    cls_label = f"cls_{unit}"
    avg_cell_area = f"avg_cell_area_{unit_area_label}"

    df = pd.DataFrame(
        {
            "resolution": resolutions,
            "number_of_cells": num_cells_list,
            avg_edge_len: avg_edge_lens,
            avg_cell_area: avg_cell_areas,
            cls_label: cls_list,
        },
        index=None,
    )

    return df


def a5stats_cli():
    """
    Command-line interface for generating A5 DGGS statistics.

    CLI options:
      -unit, --unit {m,km}
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()
    unit = args.unit
    # Get the DataFrame
    df = a5stats(unit=unit)
    # Display the DataFrame
    print(df)


def a5inspect(
    resolution: int, options={"segments": 100}, split_antimeridian: bool = False
):
    """
    Generate comprehensive inspection data for A5 DGGS cells at a given resolution.

    This function creates a detailed analysis of A5 cells including area variations,
    compactness measures, and Antimeridian crossing detection.

    Args:
        resolution: A5 resolution level (0-29)
        options: Optional dictionary of options for grid generation
        split_antimeridian: When True, apply antimeridian splitting to the resulting polygons.
            Defaults to False when None or omitted.

    Returns:
        geopandas.GeoDataFrame: DataFrame containing A5 cell inspection data with columns:
            - a5: A5 cell ID
            - resolution: Resolution level
            - geometry: Cell geometry
            - cell_area: Cell area in square meters
            - cell_perimeter: Cell perimeter in meters
            - crossed: Whether cell crosses the Antimeridian
            - norm_area: Normalized area (cell_area / mean_area)
            - ipq: Isoperimetric Quotient compactness
            - zsc: Zonal Standardized Compactness
    """
    a5_gdf = a5grid(
        resolution,
        output_format="gpd",
        options=options,
        split_antimeridian=split_antimeridian,
    )
    a5_gdf["crossed"] = a5_gdf["geometry"].apply(check_crossing_geom)
    a5_gdf = a5_gdf[~a5_gdf["crossed"]]  # remove cells that cross the Antimeridian

    mean_area = a5_gdf["cell_area"].mean()
    # Calculate normalized area
    a5_gdf["norm_area"] = a5_gdf["cell_area"] / mean_area
    # Calculate IPQ compactness using the standard formula: CI = 4πA/P²
    a5_gdf["ipq"] = 4 * np.pi * a5_gdf["cell_area"] / (a5_gdf["cell_perimeter"] ** 2)
    # Calculate zonal standardized compactness
    a5_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * a5_gdf["cell_area"]
            - np.power(a5_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / a5_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = a5_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    a5_gdf_lambert = get_cells_area(a5_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    a5_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        a5_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    a5_gdf["cvh"] = a5_gdf["cvh"].replace([np.inf, -np.inf], np.nan)
    return a5_gdf


def a5_norm_area(a5_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot normalized area map for A5 cells.

    This function creates a visualization showing how A5 cell areas vary relative
    to the mean area across the globe, highlighting areas of distortion.

    Args:
        a5_gdf: GeoDataFrame from a5inspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = a5_gdf["norm_area"].min(), 1.0, a5_gdf["norm_area"].max()
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # a5_gdf = a5_gdf[~a5_gdf["crossed"]]  # remove cells that cross the Antimeridian
    a5_gdf.to_crs(crs).plot(
        column="norm_area",
        ax=ax,
        norm=norm,
        legend=True,
        cax=cax,
        cmap="RdYlBu_r",
        legend_kwds={"label": "cell area/mean cell area", "orientation": "horizontal"},
    )
    world_countries = gpd.read_file(
        "https://raw.githubusercontent.com/opengeoshub/vopendata/refs/heads/main/shape/world_countries.geojson"
    )
    world_countries.boundary.to_crs(crs).plot(
        color=None, edgecolor="black", linewidth=0.2, ax=ax
    )
    ax.axis("off")
    cb_ax = fig.axes[1]
    cb_ax.tick_params(labelsize=14)
    cb_ax.set_xlabel(xlabel="A5 Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def a5_norm_area_hist(a5_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for A5 cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for A5 cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
        a5_gdf: GeoDataFrame from a5inspect function
    """
    # Filter out cells that cross the Antimeridian
    # a5_gdf = a5_gdf[~a5_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        a5_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (
        a5_gdf["norm_area"].min(),
        1.0,
        a5_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    # Apply colors to histogram bars using the same color mapping as the map
    for i, patch in enumerate(patches):
        # Use the center of each bin for color mapping
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.RdYlBu_r(norm(bin_center))
        patch.set_facecolor(color)

    # Add reference line at mean area (norm_area = 1)
    ax.axvline(
        x=1, color="red", linestyle="--", linewidth=2, label="Mean Area (norm_area = 1)"
    )

    # Add statistics text box
    stats_text = f"Mean: {a5_gdf['norm_area'].mean():.6f}\nStd: {a5_gdf['norm_area'].std():.6f}\nMin: {a5_gdf['norm_area'].min():.6f}\nMax: {a5_gdf['norm_area'].max():.6f}"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Customize the plot
    ax.set_xlabel("A5 normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def a5_compactness_ipq(a5_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot IPQ compactness map for A5 cells.

    This function creates a visualization showing the Isoperimetric Quotient (IPQ)
    compactness of A5 cells across the globe. IPQ measures how close each cell
    is to being circular, with values closer to 0.907 indicating more regular hexagons.

    Args:
        a5_gdf: GeoDataFrame from a5inspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # vmin, vmax, vcenter = a5_gdf['ipq'].min(), a5_gdf['ipq'].max(),np.mean([a5_gdf['ipq'].min(), a5_gdf['ipq'].max()])
    norm = TwoSlopeNorm(vmin=VMIN_PEN, vcenter=VCENTER_PEN, vmax=VMAX_PEN)
    # a5_gdf = a5_gdf[~a5_gdf["crossed"]]  # remove cells that cross the Antimeridian

    a5_gdf.to_crs(crs).plot(
        column="ipq",
        ax=ax,
        norm=norm,
        legend=True,
        cax=cax,
        cmap="viridis",
        legend_kwds={"orientation": "horizontal"},
    )
    world_countries = gpd.read_file(
        "https://raw.githubusercontent.com/opengeoshub/vopendata/refs/heads/main/shape/world_countries.geojson"
    )
    world_countries.boundary.to_crs(crs).plot(
        color=None, edgecolor="black", linewidth=0.2, ax=ax
    )
    ax.axis("off")
    cb_ax = fig.axes[1]
    cb_ax.tick_params(labelsize=14)
    cb_ax.set_xlabel(xlabel="A5 IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def a5_compactness_ipq_hist(a5_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for A5 cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for A5 cells, helping
    to understand how close cells are to being regular hexagons.

    Args:
        a5_gdf: GeoDataFrame from a5inspect function
    """
    # Filter out cells that cross the Antimeridian
    # a5_gdf = a5_gdf[~a5_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        a5_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    norm = TwoSlopeNorm(vmin=VMIN_PEN, vcenter=VCENTER_PEN, vmax=VMAX_PEN)

    # Apply colors to histogram bars using the same color mapping as the map
    for i, patch in enumerate(patches):
        # Use the center of each bin for color mapping
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    # Add reference line at ideal hexagon IPQ value (0.907)
    ax.axvline(
        x=0.907,
        color="red",
        linestyle="--",
        linewidth=2,
        label="Ideal Hexagon (IPQ = 0.907)",
    )

    # Add statistics text box
    stats_text = f"Mean: {a5_gdf['ipq'].mean():.6f}\nStd: {a5_gdf['ipq'].std():.6f}\nMin: {a5_gdf['ipq'].min():.6f}\nMax: {a5_gdf['ipq'].max():.6f}"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Customize the plot
    ax.set_xlabel("A5 IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def a5_compactness_cvh(a5_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot CVH (cell area / convex hull area) compactness map for A5 cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # a5_gdf = a5_gdf[~a5_gdf["crossed"]]  # remove cells that cross the Antimeridian
    a5_gdf = a5_gdf[np.isfinite(a5_gdf["cvh"])]
    a5_gdf = a5_gdf[a5_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    a5_gdf.to_crs(crs).plot(
        column="cvh",
        ax=ax,
        norm=norm,
        legend=True,
        cax=cax,
        cmap="viridis",
        legend_kwds={"orientation": "horizontal"},
    )
    world_countries = gpd.read_file(
        "https://raw.githubusercontent.com/opengeoshub/vopendata/refs/heads/main/shape/world_countries.geojson"
    )
    world_countries.boundary.to_crs(crs).plot(
        color=None, edgecolor="black", linewidth=0.2, ax=ax
    )
    ax.axis("off")
    cb_ax = fig.axes[1]
    cb_ax.tick_params(labelsize=14)
    cb_ax.set_xlabel(xlabel="A5 CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def a5_compactness_cvh_hist(a5_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for A5 cells.
    """
    # Filter out cells that cross the Antimeridian
    #  a5_gdf = a5_gdf[~a5_gdf["crossed"]]
    a5_gdf = a5_gdf[np.isfinite(a5_gdf["cvh"])]
    a5_gdf = a5_gdf[a5_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        a5_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Color mapping centered at 1
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    # Reference line at ideal compactness
    ax.axvline(x=1, color="red", linestyle="--", linewidth=2, label="Ideal (cvh = 1)")

    stats_text = (
        f"Mean: {a5_gdf['cvh'].mean():.6f}\n"
        f"Std: {a5_gdf['cvh'].std():.6f}\n"
        f"Min: {a5_gdf['cvh'].min():.6f}\n"
        f"Max: {a5_gdf['cvh'].max():.6f}"
    )
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    ax.set_xlabel("A5 CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def a5inspect_cli():
    """
    Command-line interface for A5 cell inspection.

    CLI options:
      -r, --resolution: A5 resolution level (0-29)
      -split, --split_antimeridian: Enable antimeridian splitting (default: enabled)
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,  # default is False to avoid splitting the Antimeridian by default
        help="Enable antimeridian splitting",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
        "Example: '{\"segments\": 1000}'",
    )
    args = parser.parse_args()
    resolution = args.resolution

    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    print(
        a5inspect(
            resolution, options=options, split_antimeridian=args.split_antimeridian
        )
    )


if __name__ == "__main__":
    a5stats_cli()

"""
This module provides functions for generating statistics for ISEA3H DGGS cells.
"""

import math
import pandas as pd
import numpy as np
import argparse
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm
from vgrid.utils.constants import (
    AUTHALIC_AREA,
    DGGS_TYPES,
    VMIN_HEX,
    VMAX_HEX,
    VCENTER_HEX,
)
from vgrid.generator.isea3hgrid import isea3hgrid
from vgrid.utils.io import add_verbose_argument
from vgrid.utils.geometry import (
    check_crossing_geom,
    characteristic_length_scale,
    convexhull_from_lambert,
    get_area_perimeter_from_lambert,
    get_cells_area,
)

min_res = DGGS_TYPES["isea3h"]["min_res"]
max_res = DGGS_TYPES["isea3h"]["max_res"]


def isea3h_metrics(resolution, unit: str = "m"):  # length unit is km, area unit is km2
    """
    Calculate metrics for ISEA3H DGGS cells at a given resolution.

    Args:
        resolution: Resolution level (0-40)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, edge_length_in_unit, cell_area_in_unit_squared)
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    num_cells = 10 * (3**resolution) + 2
    avg_cell_area = AUTHALIC_AREA / num_cells  # cell area in km²
    avg_edge_len = math.sqrt(
        (2 * avg_cell_area) / (3 * math.sqrt(3))
    )  # edge length in km
    cls = characteristic_length_scale(avg_cell_area, unit=unit)

    if resolution == 0:  # icosahedron faces
        avg_edge_len = math.sqrt((4 * avg_cell_area) / math.sqrt(3))

    # Convert to requested unit
    if unit == "km":
        avg_edge_len = avg_edge_len / (10**3)
        avg_cell_area = avg_cell_area / (10**6)

    return num_cells, avg_edge_len, avg_cell_area, cls


def isea3hstats(unit: str = "m"):  # length unit is km, area unit is km2
    """
    Generate statistics for ISEA3H DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing ISEA3H DGGS statistics with columns:
            - Resolution: Resolution level (0-40)
            - Number_of_Cells: Number of cells at each resolution
            - Avg_Edge_Length_{unit}: Average edge length in the given unit
            - Avg_Cell_Area_{unit}2: Average cell area in the squared unit
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
    for res in range(min_res, max_res + 1):
        num_cells, avg_edge_len, avg_cell_area, cls = isea3h_metrics(
            res, unit=unit
        )  # length unit is km, area unit is km2
        resolutions.append(res)
        num_cells_list.append(num_cells)
        avg_edge_lens.append(avg_edge_len)
        avg_cell_areas.append(avg_cell_area)
        cls_list.append(cls)
    # Create DataFrame
    # Build column labels with unit awareness
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


def isea3hstats_cli():
    """
    Command-line interface for generating ISEA3H DGGS statistics.

    CLI options:
      -unit, --unit {m,km}
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()  # type: ignore

    unit = args.unit

    # Get the DataFrame
    df = isea3hstats(unit=unit)

    # Display the DataFrame
    print(df)


def isea3hinspect(resolution: int, fix_antimeridian: None = None, verbose=True):
    """
    Generate comprehensive inspection data for ISEA3H DGGS cells at a given resolution.

    This function creates a detailed analysis of ISEA3H cells including area variations,
    compactness measures, and Antimeridian crossing detection.

    Args:
        resolution: ISEA3H resolution level (0-40)
        fix_antimeridian: Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        verbose: Show progress bars. Defaults to True.
    Returns:
        geopandas.GeoDataFrame: DataFrame containing ISEA3H cell inspection data with columns:
            - isea3h: ISEA3H cell ID
            - resolution: Resolution level
            - geometry: Cell geometry
            - cell_area: Cell area in square meters
            - cell_perimeter: Cell perimeter in meters
            - crossed: Whether cell crosses the Antimeridian
            - norm_area: Normalized area (cell_area / mean_area)
            - ipq: Isoperimetric Quotient compactness
            - zsc: Zonal Standardized Compactness
    """
    # Allow running on all platforms

    isea3h_gdf = isea3hgrid(
        resolution, output_format="gpd", fix_antimeridian=fix_antimeridian, verbose=verbose
    )  # remove cells that cross the Antimeridian
    isea3h_gdf["crossed"] = isea3h_gdf["geometry"].apply(check_crossing_geom)
    isea3h_gdf = isea3h_gdf[
        ~isea3h_gdf["crossed"]
    ]  # remove cells that cross the Antimeridian
    mean_area = isea3h_gdf["cell_area"].mean()
    # Calculate normalized area
    isea3h_gdf["norm_area"] = isea3h_gdf["cell_area"] / mean_area
    # Calculate IPQ compactness using the standard formula: CI = 4πA/P²
    isea3h_gdf["ipq"] = (
        4 * np.pi * isea3h_gdf["cell_area"] / (isea3h_gdf["cell_perimeter"] ** 2)
    )
    # Calculate zonal standardized compactness
    isea3h_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * isea3h_gdf["cell_area"]
            - np.power(isea3h_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / isea3h_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = isea3h_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    isea3h_gdf_lambert = get_cells_area(isea3h_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    isea3h_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        isea3h_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    isea3h_gdf["cvh"] = isea3h_gdf["cvh"].replace([np.inf, -np.inf], np.nan)
    return isea3h_gdf


def isea3h_norm_area(isea3h_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot normalized area map for ISEA3H cells.

    This function creates a visualization showing how ISEA3H cell areas vary relative
    to the mean area across the globe, highlighting areas of distortion.

    Args:
        isea3h_gdf: GeoDataFrame from isea3hinspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        isea3h_gdf["norm_area"].min(),
        1.0,
        isea3h_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # isea3h_gdf = isea3h_gdf[~isea3h_gdf["crossed"]]  # remove cells that cross the Antimeridian if split_antimeridian is False
    isea3h_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="ISEA3H Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def isea3h_compactness_ipq(isea3h_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot IPQ compactness map for ISEA3H cells.

    This function creates a visualization showing the Isoperimetric Quotient (IPQ)
    compactness of ISEA3H cells across the globe. IPQ measures how close each cell
    is to being circular, with values closer to 0.907 indicating more regular hexagons.

    Args:
        isea3h_gdf: GeoDataFrame from isea3hinspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # vmin, vcenter, vmax = isea3h_gdf['ipq'].min(), isea3h_gdf['ipq'].max(), np.mean([isea3h_gdf['ipq'].min(), isea3h_gdf['ipq'].max()])
    norm = TwoSlopeNorm(vmin=VMIN_HEX, vcenter=VCENTER_HEX, vmax=VMAX_HEX)
    # isea3h_gdf = isea3h_gdf[~isea3h_gdf["crossed"]]  # remove cells that cross the Antimeridian if split_antimeridian is False
    isea3h_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="ISEA3H IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def isea3h_norm_area_hist(isea3h_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for ISEA3H cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for ISEA3H cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
        isea3h_gdf: GeoDataFrame from isea3hinspect function
    """
    # Filter out cells that cross the Antimeridian
    # isea3h_gdf = isea3h_gdf[~isea3h_gdf["crossed"]]  # remove cells that cross the Antimeridian if split_antimeridian is False

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        isea3h_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (
        isea3h_gdf["norm_area"].min(),
        1.0,
        isea3h_gdf["norm_area"].max(),
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
    stats_text = f"Mean: {isea3h_gdf['norm_area'].mean():.3f}\nStd: {isea3h_gdf['norm_area'].std():.3f}\nMin: {isea3h_gdf['norm_area'].min():.3f}\nMax: {isea3h_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel("ISEA3H normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def isea3h_compactness_ipq_hist(isea3h_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for ISEA3H cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for ISEA3H cells, helping
    to understand how close cells are to being regular hexagons.

    Args:
        isea3h_gdf: GeoDataFrame from isea3hinspect function
    """
    # Filter out cells that cross the Antimeridian
    # isea3h_gdf = isea3h_gdf[~isea3h_gdf["crossed"]]  # remove cells that cross the Antimeridian if split_antimeridian is False

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        isea3h_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    norm = TwoSlopeNorm(vmin=VMIN_HEX, vcenter=VCENTER_HEX, vmax=VMAX_HEX)

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
    stats_text = f"Mean: {isea3h_gdf['ipq'].mean():.3f}\nStd: {isea3h_gdf['ipq'].std():.3f}\nMin: {isea3h_gdf['ipq'].min():.3f}\nMax: {isea3h_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel("ISEA3H IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def isea3h_compactness_cvh(isea3h_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot CVH (cell area / convex hull area) compactness map for ISEA3H cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # isea3h_gdf = isea3h_gdf[~isea3h_gdf["crossed"]]  # remove cells that cross the Antimeridian if split_antimeridian is False
    isea3h_gdf = isea3h_gdf[np.isfinite(isea3h_gdf["cvh"])]
    isea3h_gdf = isea3h_gdf[isea3h_gdf["cvh"] <= 1.1]

    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    isea3h_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="ISEA3H CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def isea3h_compactness_cvh_hist(isea3h_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for ISEA3H cells.
    """
    # Filter out cells that cross the Antimeridian
    # isea3h_gdf = isea3h_gdf[~isea3h_gdf["crossed"]]  # remove cells that cross the Antimeridian if split_antimeridian is False
    isea3h_gdf = isea3h_gdf[np.isfinite(isea3h_gdf["cvh"])]
    isea3h_gdf = isea3h_gdf[isea3h_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        isea3h_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
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
        f"Mean: {isea3h_gdf['cvh'].mean():.6f}\n"
        f"Std: {isea3h_gdf['cvh'].std():.6f}\n"
        f"Min: {isea3h_gdf['cvh'].min():.6f}\n"
        f"Max: {isea3h_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel("ISEA3H CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def isea3hinspect_cli():
    """
    Command-line interface for ISEA3H cell inspection.

    CLI options:
      -r, --resolution: ISEA3H resolution level (0-40)
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
    add_verbose_argument(parser)
    args = parser.parse_args()  # type: ignore
    resolution = args.resolution
    split_antimeridian = args.split_antimeridian
    print(isea3hinspect(resolution, split_antimeridian=split_antimeridian, verbose=args.verbose))


if __name__ == "__main__":
    isea3hstats_cli()

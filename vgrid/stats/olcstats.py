"""
This module provides functions for generating statistics for OLC DGGS cells.
"""

import math
import pandas as pd
import numpy as np
import argparse
import geopandas as gpd
from vgrid.utils.constants import AUTHALIC_AREA, VMIN_QUAD, VMAX_QUAD, VCENTER_QUAD
from vgrid.generator.olcgrid import olcgrid
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
from vgrid.utils.io import olc_resolutions


def olc_metrics(
    resolution: int, unit: str = "m"
):  # length unit is km, area unit is km2
    """
    Calculate metrics for OLC DGGS cells.

    Args:
        resolution: Resolution level (0-15)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, avg_edge_len_in_unit, avg_cell_area_in_unit_squared)
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # Length 2 starts with 162 cells globally
    if resolution <= 10:
        num_cells = 162 * (400 ** ((resolution // 2) - 1))
    else:
        # Length > 10: start from length 10 count, multiply by 20 per extra char
        base = 162 * (400 ** ((10 // 2) - 1))  # N(10)
        extra = resolution - 10
        num_cells = base * (20**extra)

    # Calculate area in km² first
    avg_cell_area = AUTHALIC_AREA / num_cells  # area in m2
    avg_edge_len = math.sqrt(avg_cell_area)
    cls = characteristic_length_scale(avg_cell_area, unit=unit)
    # Convert to requested unit
    if unit == "km":
        avg_cell_area = avg_cell_area / (10**6)  # Convert km² to m²
        avg_edge_len = avg_edge_len / (10**3)  # Convert km to m

    return num_cells, avg_edge_len, avg_cell_area, cls


def olcstats(unit: str = "m"):  # length unit is km, area unit is km2
    """
    Generate statistics for OLC DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing OLC DGGS statistics with columns:
            - resolution: Resolution level (2,4,6,8,10,11,12,13,14,15)
            - number_of_cells: Number of cells at each resolution
            - avg_edge_len_{unit}: Average edge length in the given unit
            - avg_cell_area_{unit}2: Average cell area in the squared unit
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # Only specific resolutions are supported

    # Initialize lists to store data
    resolutions = []
    num_cells_list = []
    avg_edge_lens = []
    avg_cell_areas = []
    cls_list = []
    for res in olc_resolutions:
        num_cells, avg_edge_len, avg_cell_area, cls = olc_metrics(res, unit=unit)
        resolutions.append(res)
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


def olcstats_cli():
    """
    Command-line interface for generating OLC DGGS statistics.

    CLI options:
      -unit, --unit {m,km}
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args, _ = parser.parse_known_args()  # type: ignore

    unit = args.unit

    # Get the DataFrame
    df = olcstats(unit=unit)

    # Display the DataFrame
    print(df)


def olcinspect(resolution: int):
    """
    Generate comprehensive inspection data for OLC DGGS cells at a given resolution.

    This function creates a detailed analysis of OLC cells including area variations,
    compactness measures, and dateline crossing detection.

    Args:
        resolution: OLC resolution level (2-15)

    Returns:
        geopandas.GeoDataFrame: DataFrame containing OLC cell inspection data with columns:
            - olc: OLC cell ID
            - resolution: Resolution level
            - geometry: Cell geometry
            - cell_area: Cell area in square meters
            - cell_perimeter: Cell perimeter in meters
            - crossed: Whether cell crosses the dateline
            - norm_area: Normalized area (cell_area / mean_area)
            - ipq: Isoperimetric Quotient compactness
            - zsc: Zonal Standardized Compactness
    """
    olc_gdf = olcgrid(resolution, output_format="gpd")
    olc_gdf["crossed"] = olc_gdf["geometry"].apply(check_crossing_geom)
    mean_area = olc_gdf["cell_area"].mean()
    # Calculate normalized area
    olc_gdf["norm_area"] = olc_gdf["cell_area"] / mean_area
    # Calculate IPQ compactness using the standard formula: CI = 4πA/P²
    olc_gdf["ipq"] = 4 * np.pi * olc_gdf["cell_area"] / (olc_gdf["cell_perimeter"] ** 2)
    # Calculate zonal standardized compactness
    olc_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * olc_gdf["cell_area"]
            - np.power(olc_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / olc_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = olc_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    olc_gdf_lambert = get_cells_area(olc_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    olc_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        olc_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    olc_gdf["cvh"] = olc_gdf["cvh"].replace([np.inf, -np.inf], np.nan)

    return olc_gdf


def olc_norm_area(olc_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot normalized area map for OLC cells.

    This function creates a visualization showing how OLC cell areas vary relative
    to the mean area across the globe, highlighting areas of distortion.

    Args:
        olc_gdf: GeoDataFrame from olcinspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vmax, vcenter = olc_gdf["norm_area"].min(), olc_gdf["norm_area"].max(), 1
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # olc_gdf = olc_gdf[~olc_gdf["crossed"]]  # remove cells that cross the dateline
    olc_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="OLC Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def olc_compactness_ipq(olc_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot IPQ compactness map for OLC cells.

    This function creates a visualization showing the Isoperimetric Quotient (IPQ)
    compactness of OLC cells across the globe. IPQ measures how close each cell
    is to being circular, with values closer to 0.785 indicating more regular squares.

    Args:
        olc_gdf: GeoDataFrame from olcinspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # vmin, vmax, vcenter = olc_gdf['ipq'].min(), olc_gdf['ipq'].max(), np.mean([olc_gdf['ipq'].min(), olc_gdf['ipq'].max()])
    norm = TwoSlopeNorm(vmin=VMIN_QUAD, vcenter=VCENTER_QUAD, vmax=VMAX_QUAD)
    # olc_gdf = olc_gdf[~olc_gdf["crossed"]]  # remove cells that cross the dateline
    olc_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="OLC IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def olc_norm_area_hist(olc_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for OLC cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for OLC cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
        olc_gdf: GeoDataFrame from olcinspect function
    """
    # Filter out cells that cross the dateline
    # olc_gdf = olc_gdf[~olc_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        olc_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vmax, vcenter = (
        olc_gdf["norm_area"].min(),
        olc_gdf["norm_area"].max(),
        1,
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
    stats_text = f"Mean: {olc_gdf['norm_area'].mean():.3f}\nStd: {olc_gdf['norm_area'].std():.3f}\nMin: {olc_gdf['norm_area'].min():.3f}\nMax: {olc_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel("OLC normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def olc_compactness_ipq_hist(olc_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for OLC cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for OLC cells, helping
    to understand how close cells are to being regular squares.

    Args:
        olc_gdf: GeoDataFrame from olcinspect function
    """
    # Filter out cells that cross the dateline
    # olc_gdf = olc_gdf[~olc_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        olc_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    norm = TwoSlopeNorm(vmin=VMIN_QUAD, vcenter=VCENTER_QUAD, vmax=VMAX_QUAD)

    # Apply colors to histogram bars using the same color mapping as the map
    for i, patch in enumerate(patches):
        # Use the center of each bin for color mapping
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    # Add reference line at ideal square IPQ value (0.785)
    ax.axvline(
        x=0.785,
        color="red",
        linestyle="--",
        linewidth=2,
        label="Ideal Square (IPQ = 0.785)",
    )

    # Add statistics text box
    stats_text = f"Mean: {olc_gdf['ipq'].mean():.3f}\nStd: {olc_gdf['ipq'].std():.3f}\nMin: {olc_gdf['ipq'].min():.3f}\nMax: {olc_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel("OLC IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def olc_compactness_cvh(olc_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot CVH (cell area / convex hull area) compactness map for ISEA4T cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # olc_gdf = olc_gdf[~olc_gdf["crossed"]]  # remove cells that cross the dateline
    olc_gdf = olc_gdf[np.isfinite(olc_gdf["cvh"])]
    olc_gdf = olc_gdf[olc_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    olc_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="OLC CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def olc_compactness_cvh_hist(olc_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for ISEA4T cells.
    """
    # Filter out cells that cross the dateline
    # olc_gdf = olc_gdf[~olc_gdf["crossed"]]
    olc_gdf = olc_gdf[np.isfinite(olc_gdf["cvh"])]
    olc_gdf = olc_gdf[olc_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        olc_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
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
        f"Mean: {olc_gdf['cvh'].mean():.6f}\n"
        f"Std: {olc_gdf['cvh'].std():.6f}\n"
        f"Min: {olc_gdf['cvh'].min():.6f}\n"
        f"Max: {olc_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel("OLC CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def olcinspect_cli():
    """
    Command-line interface for OLC cell inspection.

    CLI options:
      -r, --resolution: OLC resolution level (2-15)
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    args, _ = parser.parse_known_args()  # type: ignore
    resolution = args.resolution
    print(olcinspect(resolution))


if __name__ == "__main__":
    olcstats_cli()

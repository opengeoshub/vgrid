"""
This module provides functions for generating statistics for QTM DGGS cells.
"""

import math
import pandas as pd
import argparse
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm
from vgrid.utils.constants import (
    AUTHALIC_AREA,
    DGGS_TYPES,
    VMIN_TRI,
    VMAX_TRI,
    VCENTER_TRI,
)
from vgrid.generator.qtmgrid import qtm_grid
from vgrid.utils.geometry import (
    check_crossing_geom,
    characteristic_length_scale,
    convexhull_from_lambert,
    get_area_perimeter_from_lambert,
    get_cells_area,
)

min_res = DGGS_TYPES["qtm"]["min_res"]
max_res = DGGS_TYPES["qtm"]["max_res"]


def qtm_metrics(
    resolution: int, unit: str = "m"
):  # length unit is km, area unit is km2
    """
    Calculate metrics for QTM DGGS cells.

    Args:
        resolution: Resolution level (1-24)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, avg_edge_len_in_unit, avg_cell_area_in_unit_squared)
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    num_cells = 8 * 4 ** (resolution - 1)

    avg_cell_area = AUTHALIC_AREA / num_cells  # area in m2
    avg_edge_len = math.sqrt((4 * avg_cell_area) / math.sqrt(3))
    cls = characteristic_length_scale(avg_cell_area, unit=unit)
    # Convert to requested unit
    if unit == "km":
        avg_cell_area = avg_cell_area / (10**6)  # Convert km² to m²
        avg_edge_len = avg_edge_len / (10**3)  # Convert km to m

    return num_cells, avg_edge_len, avg_cell_area, cls


def qtmstats(unit: str = "m"):  # length unit is km, area unit is km2
    """
    Generate statistics for QTM DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing QTM DGGS statistics with columns:
            - resolution: Resolution level (1-24)
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
    for res in range(min_res, max_res + 1):
        num_cells, avg_edge_len, avg_cell_area, cls = qtm_metrics(res, unit=unit)
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


def qtmstats_cli():
    """
    Command-line interface for generating QTM DGGS statistics.

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
    df = qtmstats(unit=unit)
    # Display the DataFrame
    print(df)


def qtminspect(resolution: int):
    """
    Generate comprehensive inspection data for QTM DGGS cells at a given resolution.

    This function creates a detailed analysis of QTM cells including area variations,
    compactness measures, and dateline crossing detection.

    Args:
        resolution: QTM resolution level (1-24)

    Returns:
        geopandas.GeoDataFrame: DataFrame containing QTM cell inspection data with columns:
            - qtm: QTM cell ID
            - resolution: Resolution level
            - geometry: Cell geometry
            - cell_area: Cell area in square meters
            - cell_perimeter: Cell perimeter in meters
            - crossed: Whether cell crosses the dateline
            - norm_area: Normalized area (cell_area / mean_area)
            - ipq: Isoperimetric Quotient compactness
            - zsc: Zonal Standardized Compactness
    """
    qtm_gdf = qtm_grid(resolution)
    qtm_gdf["crossed"] = qtm_gdf["geometry"].apply(check_crossing_geom)
    qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]  # remove cells that cross the Antimeridian
    # mean_area = qtm_gdf["cell_area"].mean()
    num_cells = 8 * 4 ** (resolution - 1)
    mean_area = AUTHALIC_AREA / num_cells
    # Calculate normalized area
    qtm_gdf["norm_area"] = qtm_gdf["cell_area"] / mean_area
    # Calculate IPQ compactness using the standard formula: CI = 4πA/P²
    qtm_gdf["ipq"] = 4 * np.pi * qtm_gdf["cell_area"] / (qtm_gdf["cell_perimeter"] ** 2)
    # Calculate zonal standardized compactness
    qtm_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * qtm_gdf["cell_area"]
            - np.power(qtm_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / qtm_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = qtm_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    qtm_gdf_lambert = get_cells_area(qtm_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    qtm_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        qtm_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    qtm_gdf["cvh"] = qtm_gdf["cvh"].replace([np.inf, -np.inf], np.nan)

    return qtm_gdf


def qtm_norm_area(qtm_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot normalized area map for QTM cells.

    This function creates a visualization showing how QTM cell areas vary relative
    to the mean area across the globe, highlighting areas of distortion.

    Args:
        qtm_gdf: GeoDataFrame from qtminspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        qtm_gdf["norm_area"].min(),
        1.0,
        qtm_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]  # remove cells that cross the dateline
    qtm_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="QTM Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def qtm_norm_area_hist(qtm_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for QTM cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for QTM cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
        qtm_gdf: GeoDataFrame from qtminspect function
    """
    # Filter out cells that cross the dateline
    # qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        qtm_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (qtm_gdf["norm_area"].min(), 1.0, qtm_gdf["norm_area"].max())
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
    stats_text = f"Mean: {qtm_gdf['norm_area'].mean():.3f}\nStd: {qtm_gdf['norm_area'].std():.3f}\nMin: {qtm_gdf['norm_area'].min():.3f}\nMax: {qtm_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel("QTM normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def qtm_compactness_ipq(qtm_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot IPQ compactness map for QTM cells.

    This function creates a visualization showing the Isoperimetric Quotient (IPQ)
    compactness of QTM cells across the globe. IPQ measures how close each cell
    is to being circular, with values closer to 0.907 indicating more regular hexagons.

    Args:
        qtm_gdf: GeoDataFrame from qtminspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # vmin, vmax, vcenter = qtm_gdf['ipq'].min(), qtm_gdf['ipq'].max(), np.mean([qtm_gdf['ipq'].min(), qtm_gdf['ipq'].max()])
    norm = TwoSlopeNorm(vmin=VMIN_TRI, vcenter=VCENTER_TRI, vmax=VMAX_TRI)
    # qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]  # remove cells that cross the dateline
    qtm_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="QTM IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def qtm_compactness_ipq_hist(qtm_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for QTM cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for QTM cells, helping
    to understand how close cells are to being regular triangles.

    Args:
        qtm_gdf: GeoDataFrame from qtminspect function
    """
    # Filter out cells that cross the dateline
    # qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        qtm_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    norm = TwoSlopeNorm(vmin=VMIN_TRI, vcenter=VCENTER_TRI, vmax=VMAX_TRI)

    # Apply colors to histogram bars using the same color mapping as the map
    for i, patch in enumerate(patches):
        # Use the center of each bin for color mapping
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    # Add reference line at ideal triangle IPQ value (0.604)
    ax.axvline(
        x=0.604,
        color="red",
        linestyle="--",
        linewidth=2,
        label="Ideal Triangle (IPQ = 0.604)",
    )

    # Add statistics text box
    stats_text = f"Mean: {qtm_gdf['ipq'].mean():.3f}\nStd: {qtm_gdf['ipq'].std():.3f}\nMin: {qtm_gdf['ipq'].min():.3f}\nMax: {qtm_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel("QTM IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def qtm_compactness_cvh(qtm_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot CVH (cell area / convex hull area) compactness map for ISEA4T cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]  # remove cells that cross the dateline
    qtm_gdf = qtm_gdf[np.isfinite(qtm_gdf["cvh"])]
    qtm_gdf = qtm_gdf[qtm_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    qtm_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="QTM CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def qtm_compactness_cvh_hist(qtm_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for ISEA4T cells.
    """
    # Filter out cells that cross the dateline
    #  qtm_gdf = qtm_gdf[~qtm_gdf["crossed"]]
    qtm_gdf = qtm_gdf[np.isfinite(qtm_gdf["cvh"])]
    qtm_gdf = qtm_gdf[qtm_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        qtm_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
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
        f"Mean: {qtm_gdf['cvh'].mean():.6f}\n"
        f"Std: {qtm_gdf['cvh'].std():.6f}\n"
        f"Min: {qtm_gdf['cvh'].min():.6f}\n"
        f"Max: {qtm_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel("QTM CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def qtminspect_cli():
    """
    Command-line interface for QTM cell inspection.

    CLI options:
      -r, --resolution: QTM resolution level (1-24)
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    args = parser.parse_args()
    resolution = args.resolution
    print(qtminspect(resolution))


if __name__ == "__main__":
    qtmstats_cli()

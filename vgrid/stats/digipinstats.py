"""
This module provides functions for generating statistics for DIGIPIN DGGS cells.
"""

import math
import pandas as pd
import numpy as np
import argparse
import geopandas as gpd
from vgrid.utils.constants import (
    DGGS_TYPES,
    VMIN_QUAD,
    VMAX_QUAD,
    VCENTER_QUAD,
    AUTHALIC_RADIUS,
)
from vgrid.generator.digipingrid import digipingrid
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
from vgrid.dggs.digipin import BOUNDS
from vgrid.utils.io import add_verbose_argument, validate_digipin_resolution

min_res = DGGS_TYPES["digipin"]["min_res"]
max_res = DGGS_TYPES["digipin"]["max_res"]

# DIGIPIN regional area in square meters (India region: ~36° x 36°)
# Approximate area of India: ~3.287 million km² = 3.287 × 10^12 m²
# For calculation purposes, we use the rectangular bounds
LAT_RANGE = BOUNDS["maxLat"] - BOUNDS["minLat"]  # 36 degrees
LON_RANGE = BOUNDS["maxLon"] - BOUNDS["minLon"]  # 36 degrees

# Calculate approximate area in m² using spherical approximation
# Area ≈ (lon_range * cos(mean_lat) * R) × (lat_range * R) where R is Earth radius
MEAN_LAT = (BOUNDS["minLat"] + BOUNDS["maxLat"]) / 2
DIGIPIN_AREA = (
    LON_RANGE * math.pi / 180 * AUTHALIC_RADIUS * math.cos(MEAN_LAT * math.pi / 180)
) * (LAT_RANGE * math.pi / 180 * AUTHALIC_RADIUS)


def digipin_metrics(resolution, unit: str = "m"):  # length unit is km, area unit is km2
    """
    Calculate metrics for DIGIPIN DGGS cells.

    Args:
        resolution: Resolution level (1-10)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, avg_edge_len_in_unit, avg_cell_area_in_unit_squared, cls)
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # DIGIPIN grid has 16 (4x4) cells at each level
    # Each subdivision adds 16 cells per parent cell
    # At resolution 1: 16 cells (4x4 grid over the region)
    num_cells = 16**resolution

    # Calculate area in m² first
    avg_cell_area = DIGIPIN_AREA / num_cells
    avg_edge_len = math.sqrt(avg_cell_area)
    cls = characteristic_length_scale(avg_cell_area, unit=unit)

    # Convert to requested unit
    if unit == "km":
        avg_cell_area = avg_cell_area / (10**6)  # Convert m² to km²
        avg_edge_len = avg_edge_len / (10**3)  # Convert m to km

    return num_cells, avg_edge_len, avg_cell_area, cls


def digipinstats(unit: str = "m"):  # length unit is km, area unit is km2
    """
    Generate statistics for DIGIPIN DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing DIGIPIN DGGS statistics with columns:
            - resolution: Resolution level (1-10)
            - number_of_cells: Number of cells at each resolution
            - avg_edge_len_{unit}: Average edge length in the given unit
            - avg_cell_area_{unit}2: Average cell area in the squared unit
            - cls_{unit}: Characteristic length scale in the given unit
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
        num_cells, avg_edge_len, avg_cell_area, cls = digipin_metrics(
            res, unit=unit
        )  # length unit is km, area unit is km2
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


def digipinstats_cli():
    """
    Command-line interface for generating DIGIPIN DGGS statistics.

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
    df = digipinstats(unit=unit)

    # Display the DataFrame
    print(df)


def digipininspect(resolution, verbose=True):
    """
    Generate comprehensive inspection data for DIGIPIN DGGS cells at a given resolution.

    This function creates a detailed analysis of DIGIPIN cells including area variations,
    compactness measures, and Antimeridian crossing detection.

    Args:
        resolution: DIGIPIN resolution level (1-10)

        verbose: Show progress bars. Defaults to True.
    Returns:
        geopandas.GeoDataFrame: DataFrame containing DIGIPIN cell inspection data with columns:
            - digipin: DIGIPIN cell ID
            - resolution: Resolution level
            - geometry: Cell geometry
            - cell_area: Cell area in square meters
            - cell_perimeter: Cell perimeter in meters
            - crossed: Whether cell crosses the Antimeridian
            - norm_area: Normalized area (cell_area / mean_area)
            - ipq: Isoperimetric Quotient compactness
            - zsc: Zonal Standardized Compactness
            - cvh: Convex Hull compactness
    """
    resolution = validate_digipin_resolution(resolution)
    digipin_gdf = digipingrid(resolution, output_format="gpd", verbose=verbose)
    digipin_gdf["crossed"] = digipin_gdf["geometry"].apply(check_crossing_geom)
    # mean_area = digipin_gdf["cell_area"].mean()
    num_cells = 16**resolution
    mean_area = DIGIPIN_AREA / num_cells
    # Calculate normalized area
    digipin_gdf["norm_area"] = digipin_gdf["cell_area"] / mean_area
    # Calculate IPQ compactness using the standard formula: CI = 4πA/P²
    digipin_gdf["ipq"] = (
        4 * np.pi * digipin_gdf["cell_area"] / (digipin_gdf["cell_perimeter"] ** 2)
    )
    # Calculate zonal standardized compactness
    digipin_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * digipin_gdf["cell_area"]
            - np.power(digipin_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / digipin_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = digipin_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    digipin_gdf_lambert = get_cells_area(digipin_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    digipin_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        digipin_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    digipin_gdf["cvh"] = digipin_gdf["cvh"].replace([np.inf, -np.inf], np.nan)

    return digipin_gdf


def digipin_norm_area(digipin_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot normalized area map for DIGIPIN cells.

    This function creates a visualization showing how DIGIPIN cell areas vary relative
    to the mean area across the India region, highlighting areas of distortion.

    Args:
        digipin_gdf: GeoDataFrame from digipininspect function
        crs: Coordinate reference system for plotting (default: Mollweide projection)
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        digipin_gdf["norm_area"].min(),
        1.0,
        digipin_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # digipin_gdf = digipin_gdf[~digipin_gdf["crossed"]]  # remove cells that cross the Antimeridian
    digipin_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="DIGIPIN Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def digipin_norm_area_hist(digipin_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for DIGIPIN cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for DIGIPIN cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
        digipin_gdf: GeoDataFrame from digipininspect function
    """
    # Filter out cells that cross the Antimeridian
    # digipin_gdf = digipin_gdf[~digipin_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        digipin_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (
        digipin_gdf["norm_area"].min(),
        1.0,
        digipin_gdf["norm_area"].max(),
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
    stats_text = f"Mean: {digipin_gdf['norm_area'].mean():.3f}\nStd: {digipin_gdf['norm_area'].std():.3f}\nMin: {digipin_gdf['norm_area'].min():.3f}\nMax: {digipin_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel("DIGIPIN normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def digipin_compactness_ipq(
    digipin_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"
):
    """
    Plot IPQ compactness map for DIGIPIN cells.

    This function creates a visualization showing the Isoperimetric Quotient (IPQ)
    compactness of DIGIPIN cells across the India region. IPQ measures how close each cell
    is to being circular, with values closer to 0.785 indicating more regular squares.

    Args:
        digipin_gdf: GeoDataFrame from digipininspect function
        crs: Coordinate reference system for plotting (default: Mollweide projection)
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # vmin, vmax, vcenter = digipin_gdf['ipq'].min(), digipin_gdf['ipq'].max(), np.mean([digipin_gdf['ipq'].min(), digipin_gdf['ipq'].max()])
    norm = TwoSlopeNorm(vmin=VMIN_QUAD, vcenter=VCENTER_QUAD, vmax=VMAX_QUAD)
    # digipin_gdf = digipin_gdf[~digipin_gdf["crossed"]]  # remove cells that cross the Antimeridian
    digipin_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="DIGIPIN IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def digipin_compactness_ipq_hist(digipin_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for DIGIPIN cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for DIGIPIN cells, helping
    to understand how close cells are to being regular squares.

    Args:
        digipin_gdf: GeoDataFrame from digipininspect function
    """
    # Filter out cells that cross the Antimeridian
    # digipin_gdf = digipin_gdf[~digipin_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        digipin_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
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
    stats_text = f"Mean: {digipin_gdf['ipq'].mean():.3f}\nStd: {digipin_gdf['ipq'].std():.3f}\nMin: {digipin_gdf['ipq'].min():.3f}\nMax: {digipin_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel("DIGIPIN IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def digipin_compactness_cvh(
    digipin_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"
):
    """
    Plot CVH (cell area / convex hull area) compactness map for DIGIPIN cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.

    Args:
        digipin_gdf: GeoDataFrame from digipininspect function
        crs: Coordinate reference system for plotting (default: Mollweide projection)
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # digipin_gdf = digipin_gdf[~digipin_gdf["crossed"]]  # remove cells that cross the Antimeridian
    digipin_gdf = digipin_gdf[np.isfinite(digipin_gdf["cvh"])]
    digipin_gdf = digipin_gdf[digipin_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    digipin_gdf.to_crs(crs).plot(
        column="cvh",
        ax=ax,
        norm=norm,
        legend=True,
        cax=cax,
        cmap="viridis",
        legend_kwds={"orientation": "horizontal"},
    )
    world_countries = gpd.read_file(
        "https://raw.githubusercontent.com/opengeoshub/vopendata/refs/heads/main/shape/world_countries.geojson",
    )
    world_countries.boundary.to_crs(crs).plot(
        color=None, edgecolor="black", linewidth=0.2, ax=ax
    )
    ax.axis("off")
    cb_ax = fig.axes[1]
    cb_ax.tick_params(labelsize=14)
    cb_ax.set_xlabel(xlabel="DIGIPIN CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def digipin_compactness_cvh_hist(digipin_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for DIGIPIN cells.

    Args:
        digipin_gdf: GeoDataFrame from digipininspect function
    """
    # Filter out cells that cross the Antimeridian
    # digipin_gdf = digipin_gdf[~digipin_gdf["crossed"]]
    digipin_gdf = digipin_gdf[np.isfinite(digipin_gdf["cvh"])]
    digipin_gdf = digipin_gdf[digipin_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        digipin_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
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
        f"Mean: {digipin_gdf['cvh'].mean():.6f}\n"
        f"Std: {digipin_gdf['cvh'].std():.6f}\n"
        f"Min: {digipin_gdf['cvh'].min():.6f}\n"
        f"Max: {digipin_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel("DIGIPIN CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def digipininspect_cli():
    """
    Command-line interface for DIGIPIN cell inspection.

    CLI options:
      -r, --resolution: DIGIPIN resolution level (1-10)
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=1)
    add_verbose_argument(parser)
    args, _ = parser.parse_known_args()  # type: ignore
    resolution = args.resolution
    print(digipininspect(resolution, verbose=args.verbose))


if __name__ == "__main__":
    digipinstats_cli()

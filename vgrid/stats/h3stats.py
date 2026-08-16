"""
This module provides functions for generating statistics for H3 DGGS cells.
"""

import pandas as pd
import numpy as np

# pd.set_option('display.float_format', '{:,.3f}'.format)
import argparse
import h3
import geopandas as gpd
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
from vgrid.utils.constants import DGGS_TYPES, VMIN_HEX, VMAX_HEX, VCENTER_HEX, AUTHALIC_AREA
from vgrid.generator.h3grid import h3grid

min_res = DGGS_TYPES["h3"]["min_res"]
max_res = DGGS_TYPES["h3"]["max_res"]


def h3_metrics(resolution: int, unit: str = "m"):
    """
    Return comprehensive metrics for a resolution including number of cells,
    average edge length, average area, and area extrema analysis.

    Args:
        resolution: H3 resolution (0-15)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        dict: Dictionary containing all metrics for the resolution
    """
    length_unit = unit
    if length_unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    area_unit = {"m": "m^2", "km": "km^2"}[length_unit]

    # Basic metrics
    num_cells = h3.get_num_cells(resolution)
    avg_edge_len = h3.average_hexagon_edge_length(resolution, unit=length_unit)
    avg_area = h3.average_hexagon_area(resolution, area_unit)
    # Compute CLS (Characteristic Length Scale) always in meters first
    if length_unit == "m":
        cls = characteristic_length_scale(avg_area, unit=length_unit)
    elif length_unit == "km":
        avg_area_m2 = avg_area * (10**6)
        cls = characteristic_length_scale(avg_area_m2, unit=length_unit)

    # Return CLS in requested length unit
    # Area extrema analysis
    # Precompute base (resolution 0) hex cells (exclude pentagons)
    base_hex_cells = [idx for idx in h3.get_res0_cells() if not h3.is_pentagon(idx)]

    pentagons = list(h3.get_pentagons(resolution))

    # All hex neighbors of pentagons (exclude the pentagon cell itself)
    pentagon_neighbors = []
    for p in pentagons:
        neighbors = [n for n in h3.grid_disk(p, 1) if n != p]
        pentagon_neighbors.extend(neighbors)

    # Compute areas
    # Smallest hex area among pentagon neighbors
    min_hex_area = min(
        (h3.cell_area(idx, unit=area_unit) for idx in pentagon_neighbors),
        default=float("nan"),
    )

    # Largest hex area among center children of base hex cells
    center_children = [
        idx if resolution == 0 else h3.cell_to_center_child(idx, resolution)
        for idx in base_hex_cells
    ]
    max_hex_area = max(
        (h3.cell_area(idx, unit=area_unit) for idx in center_children),
        default=float("nan"),
    )

    # Smallest pentagon area
    min_pent_area = min(
        (h3.cell_area(idx, unit=area_unit) for idx in pentagons), default=float("nan")
    )

    # Ratios
    # hex_ratio = (
    #     (max_hex_area / min_hex_area)
    #     if (min_hex_area not in (0.0, float("nan")))
    #     else float("nan")
    # )
    hex_pent_ratio = (
        (max_hex_area / min_pent_area)
        if (min_pent_area not in (0.0, float("nan")))
        else float("nan")
    )

    return {
        "resolution": resolution,
        "number_of_cells": num_cells,
        "avg_edge_len": avg_edge_len,
        "avg_area": avg_area,
        "min_area": min_pent_area,
        "max_area": max_hex_area,
        "max_min_ratio": hex_pent_ratio,
        "cls": cls,
    }


def h3stats(unit: str = "m"):
    """
    Generate comprehensive statistics for H3 DGGS cells.

    This function combines basic H3 statistics (number of cells, edge lengths, areas)
    with area extrema analysis (min/max areas and ratios).

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: DataFrame containing comprehensive H3 DGGS statistics with columns:
            - resolution: Resolution level (0-15)
            - number_of_cells: Number of cells at each resolution
            - avg_edge_len_{unit}: Average edge length in the given unit
            - avg_area_{unit}2: Average cell area in the squared unit
            - min_area_{unit}2: Minimum pentagon area
            - max_area_{unit}2: Maximum hexagon area
            - max_min_ratio: Ratio of max hexagon area to min pentagon area
    """
    # normalize and validate unit
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # Initialize lists to store data
    resolutions = []
    num_cells_list = []
    avg_edge_lens = []
    avg_areas = []
    min_areas = []
    max_areas = []
    max_min_ratios = []
    cls_list = []
    for res in range(min_res, max_res + 1):
        # Get comprehensive metrics
        metrics_data = h3_metrics(res, unit=unit)  # length unit is km, area unit is km2

        resolutions.append(res)
        num_cells_list.append(metrics_data["number_of_cells"])
        avg_edge_lens.append(metrics_data["avg_edge_len"])
        avg_areas.append(metrics_data["avg_area"])
        min_areas.append(metrics_data["min_area"])
        max_areas.append(metrics_data["max_area"])
        max_min_ratios.append(metrics_data["max_min_ratio"])
        cls_list.append(metrics_data["cls"])
    # Create DataFrame
    # Build column labels with unit awareness (lower case)
    avg_edge_len = f"avg_edge_len_{unit}"
    avg_area = f"avg_area_{unit}"
    min_area = f"min_area_{unit}"
    max_area = f"max_area_{unit}"
    cls_label = f"cls_{unit}"
    df = pd.DataFrame(
        {
            "resolution": resolutions,
            "number_of_cells": num_cells_list,
            avg_edge_len: avg_edge_lens,
            avg_area: avg_areas,
            min_area: min_areas,
            max_area: max_areas,
            "max_min_ratio": max_min_ratios,
            cls_label: cls_list,
        }
    )

    return df


def h3stats_cli():
    """
    Command-line interface for generating H3 DGGS statistics.

    CLI options:
      -unit, --unit {m,km}
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()  # type: ignore

    unit = args.unit

    df = h3stats(unit=unit)
    df["number_of_cells"] = df["number_of_cells"].apply(lambda x: "{:,.0f}".format(x))
    print(df)


def h3inspect(resolution: int, fix_antimeridian: None = None):
    """
    Generate comprehensive inspection data for H3 DGGS cells at a given resolution.

    This function creates a detailed analysis of H3 cells including area variations,
    compactness measures, and Antimeridian crossing detection.

    Args:
        resolution: H3 resolution level (0-15)
        fix_antimeridian: Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns:
        geopandas.GeoDataFrame: DataFrame containing H3 cell inspection data with columns:
            - h3: H3 cell ID
            - resolution: Resolution level
            - geometry: Cell geometry
            - cell_area: Cell area in square meters
            - cell_perimeter: Cell perimeter in meters
            - crossed: Whether cell crosses the Antimeridian
            - is_pentagon: Whether cell is a pentagon
            - norm_area: Normalized area (cell_area / mean_area)
            - ipq: Isoperimetric Quotient compactness
            - zsc: Zonal Standardized Compactness
    """
    h3_gdf = h3grid(resolution, output_format="gpd", fix_antimeridian=fix_antimeridian)
    h3_gdf["crossed"] = h3_gdf["geometry"].apply(check_crossing_geom)
    h3_gdf = h3_gdf[~h3_gdf["crossed"]]  # remove cells that cross the Antimeridian
    h3_gdf["is_pentagon"] = h3_gdf["h3"].apply(h3.is_pentagon)
    # mean_area = h3_gdf["cell_area"].mean()
    mean_area = AUTHALIC_AREA / h3.get_num_cells(resolution)
    # Calculate normalized area
    h3_gdf["norm_area"] = h3_gdf["cell_area"] / mean_area
    # Calculate IPQ compactness using the standard formula: CI = 4πA/P²
    h3_gdf["ipq"] = 4 * np.pi * h3_gdf["cell_area"] / (h3_gdf["cell_perimeter"] ** 2)
    # Calculate zonal standardized compactness
    h3_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * h3_gdf["cell_area"]
            - np.power(h3_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / h3_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = h3_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    h3_gdf_lambert = get_cells_area(h3_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    h3_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        h3_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    h3_gdf["cvh"] = h3_gdf["cvh"].replace([np.inf, -np.inf], np.nan)
    return h3_gdf


def h3_norm_area(h3_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot normalized area map for H3 cells.

    This function creates a visualization showing how H3 cell areas vary relative
    to the mean area across the globe, highlighting areas of distortion.

    Args:
        h3_gdf: GeoDataFrame from h3inspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        h3_gdf["norm_area"].min(),
        h3_gdf["norm_area"].mean(),
        h3_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # h3_gdf = h3_gdf[~h3_gdf["crossed"]]  # remove cells that cross the Antimeridian
    h3_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="H3 Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def h3_compactness_ipq(h3_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot IPQ compactness map for H3 cells.

    This function creates a visualization showing the Isoperimetric Quotient (IPQ)
    compactness of H3 cells across the globe. IPQ measures how close each cell
    is to being circular, with values closer to 0.907 indicating more regular hexagons.

    Args:
        h3_gdf: GeoDataFrame from h3inspect function
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # vmin, vmax, vcenter = h3_gdf['ipq'].min(), h3_gdf['ipq'].max(),np.mean([h3_gdf['ipq'].min(), h3_gdf['ipq'].max()])
    norm = TwoSlopeNorm(vmin=VMIN_HEX, vcenter=VCENTER_HEX, vmax=VMAX_HEX)
    # h3_gdf = h3_gdf[~h3_gdf["crossed"]]  # remove cells that cross the Antimeridian
    h3_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="H3 IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def h3_norm_area_hist(h3_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for H3 cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for H3 cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
        h3_gdf: GeoDataFrame from h3inspect function
    """
    # Filter out cells that cross the Antimeridian
    # h3_gdf = h3_gdf[~h3_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        h3_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (
        h3_gdf["norm_area"].min(),
        h3_gdf["norm_area"].mean(),
        h3_gdf["norm_area"].max(),
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
    stats_text = f"Mean: {h3_gdf['norm_area'].mean():.3f}\nStd: {h3_gdf['norm_area'].std():.3f}\nMin: {h3_gdf['norm_area'].min():.3f}\nMax: {h3_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel("H3 normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    # Set y-axis ticks to 200, 400, 600 intervals
    y_max = ax.get_ylim()[1]
    y_ticks = np.arange(0, y_max + 200, 200)
    ax.set_yticks(y_ticks)

    plt.tight_layout()


def h3_compactness_ipq_hist(h3_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for H3 cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for H3 cells, helping
    to understand how close cells are to being regular hexagons.

    Args:
        h3_gdf: GeoDataFrame from h3inspect function
    """
    # Filter out cells that cross the Antimeridian
    # h3_gdf = h3_gdf[~h3_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        h3_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
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
    stats_text = f"Mean: {h3_gdf['ipq'].mean():.3f}\nStd: {h3_gdf['ipq'].std():.3f}\nMin: {h3_gdf['ipq'].min():.3f}\nMax: {h3_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel("H3 IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def h3_compactness_cvh(h3_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """
    Plot CVH (cell area / convex hull area) compactness map for H3 cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # h3_gdf = h3_gdf[~h3_gdf["crossed"]]  # remove cells that cross the Antimeridian
    h3_gdf = h3_gdf[np.isfinite(h3_gdf["cvh"])]
    h3_gdf = h3_gdf[h3_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    h3_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="H3 CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def h3_compactness_cvh_hist(h3_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for H3 cells.
    """
    # Filter out cells that cross the Antimeridian
    #  h3_gdf = h3_gdf[~h3_gdf["crossed"]]
    h3_gdf = h3_gdf[np.isfinite(h3_gdf["cvh"])]
    h3_gdf = h3_gdf[h3_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        h3_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
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
        f"Mean: {h3_gdf['cvh'].mean():.6f}\n"
        f"Std: {h3_gdf['cvh'].std():.6f}\n"
        f"Min: {h3_gdf['cvh'].min():.6f}\n"
        f"Max: {h3_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel("H3 CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def h3inspect_cli():
    """
    Command-line interface for H3 cell inspection.

    CLI options:
      -r, --resolution: H3 resolution level (0-15)
      -split, --split_antimeridian: Apply antimeridian fixing to the resulting polygons
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    parser.add_argument(
        "-fix--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()  # type: ignore
    resolution = args.resolution
    print(h3inspect(resolution, fix_antimeridian=args.fix_antimeridian))


if __name__ == "__main__":
    h3stats_cli()

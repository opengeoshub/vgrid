"""
DGGRID Statistics Module

This module provides functions to calculate and display statistics for DGGRID
Discrete Global Grid System (DGGS) types. It supports both command-line interface
and direct function calls.

Key Functions:
- dggrid_stats: Calculate and display statistics for a given DGGRID DGGS type and resolution
- dggridinspect: Generate detailed inspection data for a given DGGRID DGGS type and resolution
- main: Command-line interface for dggrid_stats
"""

import argparse
import json
import math
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm

from vgrid.utils.geometry import (
    check_crossing_geom,
    convexhull_from_lambert,
    dggrid_num_edges,
    get_area_perimeter_from_lambert,
    get_cells_area,
)
from vgrid.utils.constants import (
    AUTHALIC_AREA,
    VMIN_QUAD,
    VMAX_QUAD,
    VCENTER_QUAD,
    VMIN_HEX,
    VMAX_HEX,
    VCENTER_HEX,
)
from vgrid.generator.dggridgen import dggridgen
from dggrid4py import dggs_types
from pyproj import Geod
from vgrid.utils.io import (
    add_verbose_argument,
    validate_dggrid_type,
    validate_dggrid_resolution,
    create_dggrid_instance,
)
from vgrid.utils.constants import DGGRID_TYPES

geod = Geod(ellps="WGS84")


def dggridstats(
    dggrid_instance, dggs_type: str, unit: str = "m"
) -> pd.DataFrame:  # length unit is km, area unit is km2
    """length unit is m, area unit is m2
    Return a DataFrame of DGGRID stats per resolution.

    unit: 'km' or 'm' for length columns; area is squared unit.
          DGGRID native output is km^2 for area and km for CLS.

    Columns include ``avg_edge_len_{unit}``, derived from mean cell area
    and topology (hex vs quad).
    """
    dggs_type = validate_dggrid_type(dggs_type)
    unit_norm = unit.strip().lower()

    if unit_norm not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    max_res = int(DGGRID_TYPES[dggs_type]["max_res"])

    dggrid_stats_table = dggrid_instance.grid_stats_table(dggs_type, max_res)
    # Characteristic Length Scale (CLS): the diameter of a spherical cap of the same area as a cell of the specified res
    if isinstance(dggrid_stats_table, pd.DataFrame):
        rename_map = {
            "Resolution": "resolution",
            "Cells": "number_of_cells",
            "Area (km^2)": "area_km2",
            "CLS (km)": "cls_km",
        }
        dggrid_stats = dggrid_stats_table.rename(columns=rename_map).copy()
    else:
        dggrid_stats = pd.DataFrame(
            dggrid_stats_table,
            columns=["resolution", "number_of_cells", "area_km2", "cls_km"],
        )

    if unit_norm == "m":
        dggrid_stats = dggrid_stats.rename(
            columns={"area_km2": "area_m2", "cls_km": "cls_m"}
        )
        dggrid_stats["area_m2"] = dggrid_stats["area_m2"] * (10**6)
        dggrid_stats["cls_m"] = dggrid_stats["cls_m"] * (10**3)

    area_col = "area_m2" if unit_norm == "m" else "area_km2"
    avg_edge_len_col = f"avg_edge_len_{unit_norm}"
    areas = dggrid_stats[area_col]
    if dggrid_num_edges(dggs_type) == 6:
        dggrid_stats[avg_edge_len_col] = np.sqrt((2 * areas) / (3 * math.sqrt(3)))
    else:
        dggrid_stats[avg_edge_len_col] = np.sqrt(areas)

    if unit_norm == "m":
        col_order = [
            "resolution",
            "number_of_cells",
            avg_edge_len_col,
            "area_m2",
            "cls_m",
        ]
    else:
        col_order = [
            "resolution",
            "number_of_cells",
            avg_edge_len_col,
            "area_km2",
            "cls_km",
        ]
    return dggrid_stats[col_order]


def dggrid_metrics(
    dggrid_instance,
    dggs_type: str,
    resolution: int,
    unit: str = "m",
    grid_stats: pd.DataFrame | None = None,
):
    """
    Cell metrics for a DGGRID type and resolution from ``grid_stats_table``.

    Args:
        dggrid_instance: DGGRID instance
        dggs_type: DGGRID DGGS type
        resolution: Resolution level
        unit: ``'m'`` or ``'km'`` for length; area uses squared unit
        grid_stats: Optional precomputed table from :func:`dggridstats`

    Returns:
        tuple: (num_cells, avg_edge_len, avg_cell_area, cls) in requested units
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))

    unit_norm = unit.strip().lower()
    if unit_norm not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    if grid_stats is None:
        grid_stats = dggridstats(dggrid_instance, dggs_type, unit=unit_norm)

    row = grid_stats[grid_stats["resolution"] == resolution]
    if row.empty:
        raise ValueError(f"No DGGRID stats for {dggs_type} at resolution {resolution}")

    area_col = "area_m2" if unit_norm == "m" else "area_km2"
    cls_col = "cls_m" if unit_norm == "m" else "cls_km"
    avg_edge_len_col = f"avg_edge_len_{unit_norm}"
    num_cells = int(row["number_of_cells"].iloc[0])
    avg_cell_area = float(row[area_col].iloc[0])
    cls = float(row[cls_col].iloc[0])
    if avg_edge_len_col not in row.columns:
        raise ValueError(
            f"grid_stats missing {avg_edge_len_col!r}; "
            "use a table from dggridstats() or recompute stats."
        )
    avg_edge_len = float(row[avg_edge_len_col].iloc[0])

    return num_cells, avg_edge_len, avg_cell_area, cls


def dggridstats_cli():
    """
    Command-line interface for generating DGGAL DGGS statistics.

    CLI options:
      -dggs, --dggs_type {gnosis, isea3h, isea9r, ivea3h, ivea9r, rtea3h, rtea9r, rhealpix}
      -unit, --unit {m,km}
      --minres, --maxres
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-dggs", "--dggs_type", dest="dggs_type", choices=DGGRID_TYPES.keys()
    )
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()

    dggs_type = args.dggs_type
    unit = args.unit

    dggrid_instance = create_dggrid_instance()
    result = dggridstats(dggrid_instance, dggs_type, unit)
    if result is not None:
        print(result)


def dggridinspect(
    dggrid_instance,
    dggs_type: str,
    resolution: int,
    split_antimeridian: bool = False,
    aggregate: bool = False,
    options={"densification": 30},
    verbose=True,
) -> gpd.GeoDataFrame:
    """
    Generate detailed inspection data for a DGGRID DGGS type at a given resolution.

    Args:
        dggrid_instance: DGGRID instance for grid operations
        dggs_type: DGGS type supported by DGGRID (see dggs_types)
        resolution: Resolution level
        split_antimeridian: When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False to avoid splitting the Antimeridian by default.
        aggregate: When True, aggregate the resulting polygons. Defaults to False to avoid aggregation by default.
        options (dict, optional): Options to pass to grid_cell_polygons_for_extent.
            For example: {"densification": 2} to add densification points.
            Defaults to None.
        verbose: Show progress bars. Defaults to True.
    Returns:
        geopandas.GeoDataFrame: DataFrame containing inspection data with columns:
          - name (cell identifier from DGGRID)
          - resolution
          - geometry
          - cell_area (m^2)
          - cell_perimeter (m)
          - crossed (bool)
          - norm_area (area/mean_area)
          - ipq (4πA/P²)
          - zsc (sqrt(4πA - A²/R²)/P), with R=WGS84 a
    """

    # Generate grid using dggridgen
    dggrid_gdf = dggridgen(
        dggrid_instance,
        dggs_type,
        resolution,
        output_format="gpd",
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
        verbose=verbose,
    )

    # Remove cells with null or invalid geometry
    dggrid_gdf = dggrid_gdf.dropna(subset=["geometry"])
    dggrid_gdf = dggrid_gdf[dggrid_gdf.geometry.is_valid]

    # Add dggs_type column
    dggrid_gdf["dggs_type"] = f"dggrid_{dggs_type.lower()}"

    # Rename global_id to cell_id
    if "global_id" in dggrid_gdf.columns:
        dggrid_gdf = dggrid_gdf.rename(columns={"global_id": "cell_id"})

    # Determine whether current CRS is geographic; compute metrics accordingly
    if dggrid_gdf.crs.is_geographic:
        def _geod_area_perimeter(g):
            cell_area_perimeter = geod.geometry_area_perimeter(g)
            return pd.Series(
                {
                    "cell_area": abs(cell_area_perimeter[0]),
                    "cell_perimeter": abs(cell_area_perimeter[1]),
                }
            )

        metrics = dggrid_gdf.geometry.apply(_geod_area_perimeter)
        dggrid_gdf["cell_area"] = metrics["cell_area"]
        dggrid_gdf["cell_perimeter"] = metrics["cell_perimeter"]
        dggrid_gdf["crossed"] = dggrid_gdf.geometry.apply(check_crossing_geom)
    else:
        dggrid_gdf["cell_area"] = dggrid_gdf.geometry.area
        dggrid_gdf["cell_perimeter"] = dggrid_gdf.geometry.length
        dggrid_gdf["crossed"] = False

    # Add resolution column
    dggrid_gdf["resolution"] = resolution

    dggrid_gdf = dggrid_gdf[
        ~dggrid_gdf["crossed"]
    ]  # remove cells that cross the Antimeridian
    # Calculate normalized area
    # mean_area = dggrid_gdf["cell_area"].mean()
    num_cells, _, _, _ = dggrid_metrics(dggrid_instance, dggs_type, resolution)
    mean_area = AUTHALIC_AREA / num_cells
    dggrid_gdf["norm_area"] = (
        dggrid_gdf["cell_area"] / mean_area if mean_area and mean_area != 0 else np.nan
    )

    # Calculate compactness metrics (robust formulas avoiding division by zero)
    dggrid_gdf["ipq"] = (
        4 * np.pi * dggrid_gdf["cell_area"] / (dggrid_gdf["cell_perimeter"] ** 2)
    )
    dggrid_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * dggrid_gdf["cell_area"]
            - np.power(dggrid_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / dggrid_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = dggrid_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    dggrid_gdf_lambert = get_cells_area(dggrid_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    dggrid_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        dggrid_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    dggrid_gdf["cvh"] = dggrid_gdf["cvh"].replace([np.inf, -np.inf], np.nan)

    return dggrid_gdf


def dggrid_norm_area(
    dggs_type="DGGRID",
    dggrid_gdf: gpd.GeoDataFrame = None,
    crs: str | None = "proj=moll",
):
    """
    Plot normalized area map for DGGRID cells.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        dggrid_gdf["norm_area"].min(),
        1.0,
        dggrid_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # dggrid_gdf = dggrid_gdf[~dggrid_gdf["crossed"]]  # remove cells that cross the dateline
    dggrid_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel=f"{dggs_type.upper()} Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def dggrid_compactness_ipq(
    dggs_type: str = "DGGRID",
    dggrid_gdf: gpd.GeoDataFrame = None,
    crs: str | None = "proj=moll",
):
    """
    Plot IPQ compactness map for DGGRID cells.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)

    # Determine compactness bounds based on topology
    vmin, vcenter, vmax = VMIN_QUAD, VCENTER_QUAD, VMAX_QUAD

    dggs_type_norm = str(dggs_type).strip().lower()
    if any(hex_type in dggs_type_norm for hex_type in ["3h", "4h", "7h", "43h"]):
        vmin, vcenter, vmax = VMIN_HEX, VCENTER_HEX, VMAX_HEX

    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # Only filter out antimeridian-crossed cells when plotting in EPSG:4326
    # dggrid_gdf = dggrid_gdf[~dggrid_gdf["crossed"]]  # remove cells that cross the dateline
    dggrid_gdf_plot = dggrid_gdf.to_crs(crs) if crs else dggrid_gdf
    dggrid_gdf_plot.plot(
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
    wc_plot = world_countries.boundary.to_crs(crs)
    wc_plot.plot(color=None, edgecolor="black", linewidth=0.2, ax=ax)
    ax.axis("off")
    cb_ax = fig.axes[1]
    cb_ax.tick_params(labelsize=14)
    cb_ax.set_xlabel(xlabel=f"{dggs_type.upper()} IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def dggrid_norm_area_hist(dggs_type="DGGRID", dggrid_gdf: gpd.GeoDataFrame = None):
    """
    Plot histogram of normalized area for DGGRID cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for DGGRID cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
            gdf: GeoDataFrame from dggridinspect function
            dggs_type: DGGS type name for labeling
    """
    # Filter out cells that cross the dateline
    # dggrid_gdf = dggrid_gdf[~dggrid_gdf["crossed"]]  # remove cells that cross the dateline

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        dggrid_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (
        dggrid_gdf["norm_area"].min(),
        1.0,
        dggrid_gdf["norm_area"].max(),
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
    stats_text = f"Mean: {dggrid_gdf['norm_area'].mean():.3f}\nStd: {dggrid_gdf['norm_area'].std():.3f}\nMin: {dggrid_gdf['norm_area'].min():.3f}\nMax: {dggrid_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel(f"{dggs_type.upper()} normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def dggrid_compactness_ipq_hist(
    dggs_type="DGGRID", dggrid_gdf: gpd.GeoDataFrame = None
):
    """
    Plot histogram of IPQ compactness for DGGRID cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for DGGRID cells, helping
    to understand how close cells are to being regular shapes.

    Args:
            gdf: GeoDataFrame from dggridinspect function
            dggs_type: DGGS type name for labeling and determining ideal IPQ values
    """
    # Filter out cells that cross the dateline
    # dggrid_gdf = dggrid_gdf[~dggrid_gdf["crossed"]]  # remove cells that cross the dateline

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        dggrid_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    dggs_type_norm = str(dggs_type).strip().lower()
    if any(hex_type in dggs_type_norm for hex_type in ["3h", "4h", "7h", "43h"]):
        # Hexagonal cells
        norm = TwoSlopeNorm(vmin=VMIN_HEX, vcenter=VCENTER_HEX, vmax=VMAX_HEX)
        ideal_ipq = 0.907  # Ideal hexagon
        shape_name = "Hexagon"
    else:
        # Quadrilateral cells
        norm = TwoSlopeNorm(vmin=VMIN_QUAD, vcenter=VCENTER_QUAD, vmax=VMAX_QUAD)
        ideal_ipq = 0.785  # Ideal square (π/4)
        shape_name = "Square"

    # Apply colors to histogram bars using the same color mapping as the map
    for i, patch in enumerate(patches):
        # Use the center of each bin for color mapping
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    # Add reference line at ideal IPQ value
    ax.axvline(
        x=ideal_ipq,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Ideal {shape_name} (IPQ = {ideal_ipq:.3f})",
    )

    # Add statistics text box
    stats_text = f"Mean: {dggrid_gdf['ipq'].mean():.3f}\nStd: {dggrid_gdf['ipq'].std():.3f}\nMin: {dggrid_gdf['ipq'].min():.3f}\nMax: {dggrid_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel(f"{dggs_type.upper()} IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def dggrid_compactness_cvh(
    dggs_type="DGGRID",
    dggrid_gdf: gpd.GeoDataFrame = None,
    crs: str | None = "proj=moll",
):
    """
    Plot CVH (cell area / convex hull area) compactness map for DGGRID cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # dggrid_gdf = dggrid_gdf[~dggrid_gdf["crossed"]]  # remove cells that cross the dateline
    dggrid_gdf = dggrid_gdf[np.isfinite(dggrid_gdf["cvh"])]
    dggrid_gdf = dggrid_gdf[dggrid_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    dggrid_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="DGGRID CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def dggrid_compactness_cvh_hist(
    dggs_type="DGGRID", dggrid_gdf: gpd.GeoDataFrame = None
):
    """
    Plot histogram of CVH (cell area / convex hull area) for DGGRID cells.
    """
    # Filter out cells that cross the dateline
    # dggrid_gdf = dggrid_gdf[~dggrid_gdf["crossed"]]
    dggrid_gdf = dggrid_gdf[np.isfinite(dggrid_gdf["cvh"])]
    dggrid_gdf = dggrid_gdf[dggrid_gdf["cvh"] <= 1.1]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        dggrid_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
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
        f"Mean: {dggrid_gdf['cvh'].mean():.6f}\n"
        f"Std: {dggrid_gdf['cvh'].std():.6f}\n"
        f"Min: {dggrid_gdf['cvh'].min():.6f}\n"
        f"Max: {dggrid_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel("DGGRID CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def dggridinspect_cli():
    """
    Command-line interface for DGGRID cell inspection.

    CLI options:
      -dggs, --dggs_type: DGGS type from dggs_types
      -r, --resolution: Resolution level
      --no-split_antimeridian: Disable antimeridian fixing (default: enabled)
      --no-aggregate: Disable aggregation (default: enabled)
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-dggs", "--dggs_type", dest="dggs_type", choices=dggs_types, required=True
    )
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,  # default is False to avoid splitting the Antimeridian by default
        help="Enable antimeridian fixing",
    )
    parser.add_argument(
        "-aggregate",
        "--aggregate",
        action="store_true",
        default=False,  # default is False to avoid aggregation by default
        help="Enable aggregation",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_for_extent. "
        "Example: '{\"densification\": 2}'",
    )
    add_verbose_argument(parser)
    args = parser.parse_args()
    dggrid_instance = create_dggrid_instance()
    dggs_type = args.dggs_type
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
        dggridinspect(
            dggrid_instance,
            dggs_type,
            resolution,
            split_antimeridian=args.split_antimeridian,
            aggregate=args.aggregate,
            options=options,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    dggridstats_cli()

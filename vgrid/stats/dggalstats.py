"""
This module provides lightweight wrappers for DGGAL using the external `dgg` CLI directly.

Per request, `dggalstats` simply returns the direct output from
`dgg <dggs_type> level` without computing any additional metrics.
"""

import argparse
import shutil
import subprocess
import sys
import numpy as np
import pandas as pd
import math
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import TwoSlopeNorm

from vgrid.utils.geometry import (
    check_crossing_geom,
    characteristic_length_scale,
    convexhull_from_lambert,
    get_area_perimeter_from_lambert,
    get_cells_area,
)
from vgrid.utils.constants import (
    AUTHALIC_AREA,
    DGGAL_TYPES,
    VMIN_QUAD,
    VMAX_QUAD,
    VCENTER_QUAD,
    VMIN_HEX,
    VMAX_HEX,
    VCENTER_HEX,
)
from vgrid.generator.dggalgen import dggalgen
from pyproj import Geod
from vgrid.utils.io import validate_dggal_resolution, validate_dggal_type

# Import dggal library
from dggal import *

# Initialize dggal application
app = Application(appGlobals=globals())
pydggal_setup(app)

geod = Geod(ellps="WGS84")


def dggalinfo(dggs_type: str) -> str | None:
    """
    Return the direct stdout from `dgg <dggs_type> level`.

    Args:
            dggs_type: DGGS type supported by DGGAL (see vgrid.utils.constants.dggs_type)

    Returns:
            stdout string on success; None on failure.
    """
    dggs_type = validate_dggal_type(dggs_type)

    dgg_exe = shutil.which("dgg")
    if dgg_exe is None:
        print(
            "Error: `dgg` command not found. Please ensure the `dggal` package is installed and `dgg` is on PATH.",
            file=sys.stderr,
        )
        return None

    # Use the style: `dgg <dggs_type> level`
    cmd = [dgg_exe, dggs_type, "level"]

    try:
        completed = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        stdout = completed.stdout
    except Exception as exc:
        print(f"Failed to run {' '.join(cmd)}: {exc}", file=sys.stderr)
        return None
    # Return the textual table directly for display
    return stdout


def dggalinfo_cli():
    """
    Command-line interface for generating DGGAL DGGS statistics.

    CLI options:
      -dggs, --dggs_type {gnosis, isea3h, isea9r, ivea3h, ivea9r, rtea3h, rtea9r, rhealpix}
      -unit, --unit {m,km}
      --minres, --maxres
    """
    parser = argparse.ArgumentParser(add_help=False)
    # Positional shorthand: dggalstats isea3h
    parser.add_argument("pos_dggs_type", nargs="?", choices=DGGAL_TYPES.keys())
    # Optional flag remains supported for type
    parser.add_argument(
        "-dggs", "--dggs_type", dest="dggs_type", choices=DGGAL_TYPES.keys()
    )
    args, _ = parser.parse_known_args()

    # Resolve parameters from positional or flagged inputs
    dggs_type = args.pos_dggs_type or args.dggs_type

    if dggs_type is None:
        raise SystemExit(
            "Error: dggs_type is required. Usage: dggalstats <dggs_type> or with flag -t"
        )

    result = dggalinfo(dggs_type)
    if result is not None:
        print(result)


def dggal_metrics(
    dggs_type: str, resolution: int, unit: str = "m"
):  # length unit is m, area unit is m2
    """
    Calculate metrics for DGGAL cells at a given resolution.

    Args:
            dggs_type: DGGS type supported by DGGAL (see vgrid.utils.constants.DGGAL_TYPES)
            resolution: Resolution level (0-29)
            unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
            tuple: (num_cells, edge_length_in_unit, cell_area_in_unit_squared)
    """

    dggs_type = validate_dggal_type(dggs_type)
    resolution = validate_dggal_resolution(dggs_type, int(resolution))

    unit_norm = unit.strip().lower()
    if unit_norm not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # 'gnosis','isea4r','isea9r','isea3h','isea7h','isea7h_z7',
    # 'ivea4r','ivea9r','ivea3h','ivea7h','ivea7h_z7','rtea4r','rtea9r','rtea3h','rtea7h','rtea7h_z7','healpix','rhealpix'
    num_edges = 4
    if dggs_type in [
        "isea3h",
        "isea7h",
        "isea7h_z7",
        "ivea3h",
        "ivea7h",
        "ivea7h_z7",
        "rtea3h",
        "rtea7h",
        "rtea7h_z7",
    ]:
        num_edges = 6  # Hexagonal cells

    # Calculate number of cells using the original formulas
    # Need to be rechecked
    num_cells = 1
    if dggs_type == "gnosis":
        num_cells = (16 * (4**resolution) + 8) // 3
    elif dggs_type in ["isea3h", "ivea3h", "rtea3h"]:
        num_cells = 10 * (3**resolution) + 2
    elif dggs_type in ["isea4r", "ivea4r", "rtea4r"]:
        num_cells = 10 * (4**resolution)
    elif dggs_type in [
        "isea7h",
        "isea7h_z7",
        "ivea7h",
        "ivea7h_z7",
        "rtea7h",
        "rtea7h_z7",
    ]:
        num_cells = 10 * (7**resolution) + 2
    elif dggs_type in ["isea9r", "ivea9r", "rtea9r"]:
        num_cells = 10 * (9**resolution)
    elif dggs_type in ["healpix"]:
        num_cells = 12 * (4**resolution)
    elif dggs_type in ["rhealpix"]:
        num_cells = 6 * (9**resolution)

    avg_cell_area = AUTHALIC_AREA / num_cells  # area in m2

    # Calculate average edge length based on the number of edges
    if num_edges == 6:  # Hexagonal cells
        avg_edge_len = math.sqrt((2 * avg_cell_area) / (3 * math.sqrt(3)))
    else:  # Square or other polygonal cells
        avg_edge_len = math.sqrt(avg_cell_area)

    cls = characteristic_length_scale(avg_cell_area, unit=unit)

    # Convert to requested unit
    if unit_norm == "km":
        avg_edge_len = avg_edge_len / (10**3)
        avg_cell_area = avg_cell_area / (10**6)

    return num_cells, avg_edge_len, avg_cell_area, cls


def dggalstats(
    dggs_type: str, unit: str = "m"
) -> pd.DataFrame | None:  # length unit is km, area unit is km2
    """
    Compute and return a DataFrame of DGGAL metrics per resolution for the given type.

    Args:
            dggs_type: DGGS type supported by DGGAL (see vgrid.utils.constants.DGGAL_TYPES)
            unit: 'm' or 'km' for length; area columns will reflect the squared unit

    Returns:
            pandas DataFrame with columns for resolution, number of cells, average edge length,
            and average cell area in the requested units.
    """
    dggs_type = validate_dggal_type(dggs_type)
    min_res = int(DGGAL_TYPES[dggs_type]["min_res"])
    max_res = int(DGGAL_TYPES[dggs_type]["max_res"])

    # Initialize lists to store data
    resolutions = []
    num_cells_list = []
    avg_edge_lens = []
    avg_cell_areas = []
    cls_list = []
    for res in range(min_res, max_res + 1):
        num_cells, avg_edge_len, avg_cell_area, cls = dggal_metrics(
            dggs_type, res, unit=unit
        )  # length unit is km, area unit is km2
        resolutions.append(res)
        num_cells_list.append(num_cells)
        avg_edge_lens.append(avg_edge_len)
        avg_cell_areas.append(avg_cell_area)
        cls_list.append(cls)
    # Build column labels with unit awareness
    avg_edge_len_col = f"avg_edge_len_{unit}"
    unit_area_label = {"m": "m2", "km": "km2"}[unit]
    avg_cell_area_col = f"avg_cell_area_{unit_area_label}"
    cls_label = f"cls_{unit}"
    df = pd.DataFrame(
        {
            "resolution": resolutions,
            "number_of_cells": num_cells_list,
            avg_edge_len_col: avg_edge_lens,
            avg_cell_area_col: avg_cell_areas,
            cls_label: cls_list,
        }
    )

    return df


def dggalstats_cli():
    """
    Command-line interface for generating DGGAL DGGS statistics.

    CLI options:
      -dggs, --dggs_type {gnosis, isea3h, isea9r, ivea3h, ivea9r, rtea3h, rtea9r, rhealpix}
      -unit, --unit {m,km}
      --minres, --maxres
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-dggs", "--dggs_type", dest="dggs_type", choices=DGGAL_TYPES.keys()
    )
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()

    dggs_type = args.dggs_type
    unit = args.unit

    result = dggalstats(dggs_type, unit)
    if result is not None:
        print(result)


def dggalinspect(
    dggs_type: str, resolution: int, split_antimeridian: bool = False
) -> gpd.GeoDataFrame:
    """
    Generate detailed inspection data for a DGGAL DGGS type at a given resolution.

    Args:
        dggs_type: DGGS type supported by DGGAL
        resolution: Resolution level
        split_antimeridian: When True, apply antimeridian splitting to the resulting polygons.
            Defaults to True when None or omitted.

    Returns:
        geopandas.GeoDataFrame with columns:
          - ZoneID (as provided by DGGAL output; no renaming is performed)
          - resolution
          - geometry
          - cell_area (m^2)
          - cell_perimeter (m)
          - crossed (bool)
          - norm_area (area/mean_area)
          - ipq (4πA/P²)
          - zsc (sqrt(4πA - A²/R²)/P), with R=WGS84 a
    """
    dggal_gdf = dggalgen(
        dggs_type,
        resolution,
        output_format="gpd",
        split_antimeridian=split_antimeridian,
    )

    # Determine whether current CRS is geographic; compute metrics accordingly
    if dggal_gdf.crs.is_geographic:
        def _geod_area_perimeter(g):
            cell_area_perimeter = geod.geometry_area_perimeter(g)
            return pd.Series(
                {
                    "cell_area": abs(cell_area_perimeter[0]),
                    "cell_perimeter": abs(cell_area_perimeter[1]),
                }
            )

        metrics = dggal_gdf.geometry.apply(_geod_area_perimeter)
        dggal_gdf["cell_area"] = metrics["cell_area"]
        dggal_gdf["cell_perimeter"] = metrics["cell_perimeter"]
        dggal_gdf["crossed"] = dggal_gdf.geometry.apply(check_crossing_geom)
    else:
        dggal_gdf["cell_area"] = dggal_gdf.geometry.area
        dggal_gdf["cell_perimeter"] = dggal_gdf.geometry.length
        dggal_gdf["crossed"] = False

    dggal_gdf = dggal_gdf[
        ~dggal_gdf["crossed"]
    ]  # remove cells that cross the Antimeridian
    # mean_area = dggal_gdf["cell_area"].mean()
    num_cells, _, _, _ = dggal_metrics(dggs_type, resolution)
    mean_area = AUTHALIC_AREA / num_cells

    dggal_gdf["norm_area"] = (
        dggal_gdf["cell_area"] / mean_area if mean_area and mean_area != 0 else np.nan
    )
    # Robust formulas avoiding division by zero
    dggal_gdf["ipq"] = (
        4 * np.pi * dggal_gdf["cell_area"] / (dggal_gdf["cell_perimeter"] ** 2)
    )

    dggal_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * dggal_gdf["cell_area"]
            - np.power(dggal_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / dggal_gdf["cell_perimeter"]
    )

    # Compute convex hull using Lambert Azimuthal Equal Area projection
    convex_hull = dggal_gdf["geometry"].apply(convexhull_from_lambert)
    # Calculate convex hull area using Lambert projection
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    # Calculate cell area using Lambert projection for consistent cvh calculation
    dggal_gdf_lambert = get_cells_area(dggal_gdf.copy(), "LAEA")
    # Compute CVH safely; set to NaN where convex hull area is non-positive or invalid
    dggal_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        dggal_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    # Replace any accidental inf values with NaN
    dggal_gdf["cvh"] = dggal_gdf["cvh"].replace([np.inf, -np.inf], np.nan)
    return dggal_gdf


def dggal_norm_area_hist(dggs_type: str, dggal_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of normalized area for DGGAL cells.

    This function creates a histogram visualization showing the distribution
    of normalized areas for DGGAL cells, helping to understand area variations
    and identify patterns in area distortion.

    Args:
            gdf: GeoDataFrame from dggalinspect function
            dggs_type: DGGS type name for labeling
    """
    # Filter out cells that cross the dateline
    # dggal_gdf = dggal_gdf[~dggal_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        dggal_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    vmin, vcenter, vmax = (
        dggal_gdf["norm_area"].min(),
        1.0,
        dggal_gdf["norm_area"].max(),
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
    stats_text = f"Mean: {dggal_gdf['norm_area'].mean():.6f}\nStd: {dggal_gdf['norm_area'].std():.6f}\nMin: {dggal_gdf['norm_area'].min():.6f}\nMax: {dggal_gdf['norm_area'].max():.6f}"
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


def dggal_norm_area(
    dggs_type: str, dggal_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"
):  # type: ignore
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        dggal_gdf["norm_area"].min(),
        1.0,
        dggal_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # dggal_gdf = dggal_gdf[~dggal_gdf["crossed"]]  # remove cells that cross the dateline
    dggal_gdf.to_crs(crs).plot(
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


def dggal_compactness_ipq(
    dggs_type: str,
    dggal_gdf: gpd.GeoDataFrame,
    crs: str | None = "proj=moll",  # type: ignore
):
    """
    Plot IPQ compactness map for DGGAL cells (generic visualization).
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)

    vmin, vcenter, vmax = VMIN_QUAD, VCENTER_QUAD, VMAX_QUAD

    dggs_type_norm = str(dggs_type).strip().lower()
    if dggs_type_norm in ["isea3h", "ivea3h", "rtea3h"]:
        vmin, vcenter, vmax = VMIN_HEX, VCENTER_HEX, VMAX_HEX

    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    # Only filter out antimeridian-crossed cells when plotting in EPSG:4326
    # dggal_gdf = dggal_gdf[~dggal_gdf["crossed"]]
    gdf_plot = dggal_gdf.to_crs(crs) if crs else dggal_gdf
    gdf_plot.plot(
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


def dggal_compactness_ipq_hist(dggs_type: str, dggal_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of IPQ compactness for DGGAL cells.

    This function creates a histogram visualization showing the distribution
    of Isoperimetric Quotient (IPQ) compactness values for DGGAL cells, helping
    to understand how close cells are to being regular shapes.

    Args:
            gdf: GeoDataFrame from dggalinspect function
            dggs_type: DGGS type name for labeling and determining ideal IPQ values
    """
    # Filter out cells that cross the dateline
    # dggal_gdf = dggal_gdf[~dggal_gdf["crossed"]]

    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get histogram data
    counts, bins, patches = ax.hist(
        dggal_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )

    # Create color ramp using the same normalization as the map function
    dggs_type_norm = str(dggs_type).strip().lower()
    if dggs_type_norm in ["isea3h", "ivea3h", "rtea3h"]:
        # Hexagonal cells
        norm = TwoSlopeNorm(vmin=VMIN_HEX, vcenter=VCENTER_HEX, vmax=VMAX_HEX)
        ideal_ipq = 0.907  # Ideal hexagon
        shape_name = "Hexagon"
    else:
        # Quadrilateral cells (gnosis, isea9r, ivea9r, rtea9r, rhealpix)
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
        label=f"Ideal {shape_name} (IPQ = {ideal_ipq:.6f})",
    )

    # Add statistics text box
    stats_text = f"Mean: {dggal_gdf['ipq'].mean():.6f}\nStd: {dggal_gdf['ipq'].std():.6f}\nMin: {dggal_gdf['ipq'].min():.6f}\nMax: {dggal_gdf['ipq'].max():.6f}"
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


def dggal_compactness_cvh(
    dggs_type: str, dggal_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"
):
    """
    Plot CVH (cell area / convex hull area) compactness map for DGGAL cells.

    Values are in (0, 1], with 1 indicating the most compact (convex) shape.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    # dggal_gdf = dggal_gdf[~dggal_gdf["crossed"]]  # remove cells that cross the dateline
    dggal_gdf = dggal_gdf[np.isfinite(dggal_gdf["cvh"])]
    dggal_gdf = dggal_gdf[dggal_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    dggal_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel=f"{dggs_type.upper()} CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def dggal_compactness_cvh_hist(dggs_type: str, dggal_gdf: gpd.GeoDataFrame):
    """
    Plot histogram of CVH (cell area / convex hull area) for DGGAL cells.
    """
    # Filter out cells that cross the dateline
    # dggal_gdf = dggal_gdf[~dggal_gdf["crossed"]]
    dggal_gdf = dggal_gdf[np.isfinite(dggal_gdf["cvh"])]
    dggal_gdf = dggal_gdf[dggal_gdf["cvh"] <= 1.1]
    # Create the histogram with color ramp
    fig, ax = plt.subplots(figsize=(10, 6))

    counts, bins, patches = ax.hist(
        dggal_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
    )

    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    # Reference line at ideal compactness
    ax.axvline(x=1, color="red", linestyle="--", linewidth=2, label="Ideal (cvh = 1)")

    stats_text = (
        f"Mean: {dggal_gdf['cvh'].mean():.6f}\n"
        f"Std: {dggal_gdf['cvh'].std():.6f}\n"
        f"Min: {dggal_gdf['cvh'].min():.6f}\n"
        f"Max: {dggal_gdf['cvh'].max():.6f}"
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

    ax.set_xlabel(f"{dggs_type.upper()} CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()


def dggalinspect_cli():
    """
    Command-line interface for DGGAL cell inspection.

    CLI options:
      -t, --dggs_type
      -r, --resolution
      -split, --split_antimeridian: Enable antimeridian splitting (default: enabled)
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-t", "--dggs_type", dest="dggs_type", choices=DGGAL_TYPES.keys(), required=True
    )
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,  # default is False to avoid splitting the Antimeridian by default
        help="Enable antimeridian splitting",
    )
    args = parser.parse_args()
    dggs_type = args.dggs_type
    resolution = args.resolution
    print(
        dggalinspect(dggs_type, resolution, split_antimeridian=args.split_antimeridian)
    )


if __name__ == "__main__":
    dggalstats_cli()

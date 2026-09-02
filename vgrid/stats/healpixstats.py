"""
This module provides functions for generating statistics for HEALPix DGGS cells.
"""

import math
import pandas as pd
import numpy as np
import argparse
import geopandas as gpd
from vgrid.utils.constants import (
    AUTHALIC_AREA,
    DGGS_TYPES,
    VMIN_QUAD,
    VMAX_QUAD,
    VCENTER_QUAD,
)
from vgrid.generator.healpixgrid import healpixgrid
from vgrid.utils.io import add_verbose_argument
from vgrid.dggs.healpix import nside2npix, order2nside
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

min_res = DGGS_TYPES["healpix"]["min_res"]
max_res = DGGS_TYPES["healpix"]["max_res"]


def healpix_metrics(resolution: int, unit: str = "m"):
    """
    Calculate metrics for HEALPix DGGS cells at a given resolution.

    Args:
        resolution: Order / resolution level (0-29)
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        tuple: (num_cells, edge_length_in_unit, cell_area_in_unit_squared, cls)
    """
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    num_cells = nside2npix(order2nside(resolution))
    avg_cell_area = AUTHALIC_AREA / num_cells  # m²
    avg_edge_len = math.sqrt(avg_cell_area)
    cls = characteristic_length_scale(avg_cell_area, unit=unit)
    if unit == "km":
        avg_cell_area = avg_cell_area / (10**6)
        avg_edge_len = avg_edge_len / (10**3)
        cls = cls / (10**3)
    return num_cells, avg_edge_len, avg_cell_area, cls


def healpixstats(unit: str = "m"):
    """
    Generate statistics for HEALPix DGGS cells.

    Args:
        unit: 'm' or 'km' for length; area will be 'm^2' or 'km^2'

    Returns:
        pandas.DataFrame: HEALPix DGGS statistics by resolution.
    """
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    resolutions = []
    num_cells_list = []
    avg_edge_lens = []
    avg_cell_areas = []
    cls_list = []
    for res in range(min_res, max_res + 1):
        num_cells, avg_edge_len, avg_cell_area, cls = healpix_metrics(res, unit=unit)
        resolutions.append(res)
        num_cells_list.append(num_cells)
        avg_edge_lens.append(avg_edge_len)
        avg_cell_areas.append(avg_cell_area)
        cls_list.append(cls)

    avg_edge_len = f"avg_edge_len_{unit}"
    unit_area_label = {"m": "m2", "km": "km2"}[unit]
    avg_cell_area = f"avg_cell_area_{unit_area_label}"
    cls_label = f"cls_{unit}"
    return pd.DataFrame(
        {
            "resolution": resolutions,
            "number_of_cells": num_cells_list,
            avg_edge_len: avg_edge_lens,
            avg_cell_area: avg_cell_areas,
            cls_label: cls_list,
        }
    )


def healpixstats_cli():
    """Command-line interface for generating HEALPix DGGS statistics."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-unit", "--unit", dest="unit", choices=["m", "km"], default="m"
    )
    args = parser.parse_args()
    print(healpixstats(unit=args.unit))


def healpixinspect(resolution: int = 0, fix_antimeridian: str = None, verbose=True):
    """
    Generate comprehensive inspection data for HEALPix DGGS cells at a given resolution.

    Args:
        resolution: HEALPix order (0-29)
        fix_antimeridian: Antimeridian fixing method

        verbose: Show progress bars. Defaults to True.
    Returns:
        geopandas.GeoDataFrame with area, compactness, and dateline metrics.
    """
    healpix_gdf = healpixgrid(
        resolution, output_format="gpd", fix_antimeridian=fix_antimeridian, verbose=verbose
    )
    healpix_gdf["crossed"] = healpix_gdf["geometry"].apply(check_crossing_geom)
    healpix_gdf = healpix_gdf[~healpix_gdf["crossed"]]
    # mean_area = healpix_gdf["cell_area"].mean()
    num_cells = nside2npix(order2nside(resolution))
    mean_area = AUTHALIC_AREA / num_cells
    healpix_gdf["norm_area"] = healpix_gdf["cell_area"] / mean_area
    healpix_gdf["ipq"] = (
        4 * np.pi * healpix_gdf["cell_area"] / (healpix_gdf["cell_perimeter"] ** 2)
    )
    healpix_gdf["zsc"] = (
        np.sqrt(
            4 * np.pi * healpix_gdf["cell_area"]
            - np.power(healpix_gdf["cell_area"], 2) / np.power(6378137, 2)
        )
        / healpix_gdf["cell_perimeter"]
    )
    convex_hull = healpix_gdf["geometry"].apply(convexhull_from_lambert)
    convex_hull_area = convex_hull.apply(
        lambda g: get_area_perimeter_from_lambert(g)[0] if g is not None else np.nan
    )
    healpix_gdf_lambert = get_cells_area(healpix_gdf.copy(), "LAEA")
    healpix_gdf["cvh"] = np.where(
        (convex_hull_area > 0) & np.isfinite(convex_hull_area),
        healpix_gdf_lambert["area"] / convex_hull_area,
        np.nan,
    )
    healpix_gdf["cvh"] = healpix_gdf["cvh"].replace([np.inf, -np.inf], np.nan)
    return healpix_gdf


def healpix_norm_area(healpix_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"):
    """Plot normalized area map for HEALPix cells."""
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    vmin, vcenter, vmax = (
        healpix_gdf["norm_area"].min(),
        1.0,
        healpix_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    healpix_gdf = healpix_gdf[~healpix_gdf["crossed"]]
    healpix_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="HEALPix Normalized Area", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def healpix_compactness_ipq(
    healpix_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"
):
    """Plot IPQ compactness map for HEALPix cells."""
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    norm = TwoSlopeNorm(vmin=VMIN_QUAD, vcenter=VCENTER_QUAD, vmax=VMAX_QUAD)
    healpix_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="HEALPix IPQ Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def healpix_norm_area_hist(healpix_gdf: gpd.GeoDataFrame):
    """Plot histogram of normalized area for HEALPix cells."""
    fig, ax = plt.subplots(figsize=(10, 6))
    counts, bins, patches = ax.hist(
        healpix_gdf["norm_area"], bins=50, alpha=0.7, edgecolor="black"
    )
    vmin, vcenter, vmax = (
        healpix_gdf["norm_area"].min(),
        1.0,
        healpix_gdf["norm_area"].max(),
    )
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.RdYlBu_r(norm(bin_center))
        patch.set_facecolor(color)

    ax.axvline(
        x=1, color="red", linestyle="--", linewidth=2, label="Mean Area (norm_area = 1)"
    )
    stats_text = (
        f"Mean: {healpix_gdf['norm_area'].mean():.3f}\n"
        f"Std: {healpix_gdf['norm_area'].std():.3f}\n"
        f"Min: {healpix_gdf['norm_area'].min():.3f}\n"
        f"Max: {healpix_gdf['norm_area'].max():.3f}"
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
    ax.set_xlabel("HEALPix normalized area", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()


def healpix_compactness_ipq_hist(healpix_gdf: gpd.GeoDataFrame):
    """Plot histogram of IPQ compactness for HEALPix cells."""
    fig, ax = plt.subplots(figsize=(10, 6))
    counts, bins, patches = ax.hist(
        healpix_gdf["ipq"], bins=50, alpha=0.7, edgecolor="black"
    )
    norm = TwoSlopeNorm(vmin=VMIN_QUAD, vcenter=VCENTER_QUAD, vmax=VMAX_QUAD)
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    ax.axvline(
        x=1.0,
        color="red",
        linestyle="--",
        linewidth=2,
        label="Ideal Quadrilateral (IPQ = 1.0)",
    )
    stats_text = (
        f"Mean: {healpix_gdf['ipq'].mean():.3f}\n"
        f"Std: {healpix_gdf['ipq'].std():.3f}\n"
        f"Min: {healpix_gdf['ipq'].min():.3f}\n"
        f"Max: {healpix_gdf['ipq'].max():.3f}"
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
    ax.set_xlabel("HEALPix IPQ Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()


def healpix_compactness_cvh(
    healpix_gdf: gpd.GeoDataFrame, crs: str | None = "proj=moll"
):
    """Plot CVH (cell area / convex hull area) compactness map for HEALPix cells."""
    fig, ax = plt.subplots(figsize=(10, 5))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    healpix_gdf = healpix_gdf[np.isfinite(healpix_gdf["cvh"])]
    healpix_gdf = healpix_gdf[healpix_gdf["cvh"] <= 1.1]
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    healpix_gdf.to_crs(crs).plot(
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
    cb_ax.set_xlabel(xlabel="HEALPix CVH Compactness", fontsize=14)
    ax.margins(0)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    plt.tight_layout()


def healpix_compactness_cvh_hist(healpix_gdf: gpd.GeoDataFrame):
    """Plot histogram of CVH compactness for HEALPix cells."""
    healpix_gdf = healpix_gdf[np.isfinite(healpix_gdf["cvh"])]
    healpix_gdf = healpix_gdf[healpix_gdf["cvh"] <= 1.1]

    fig, ax = plt.subplots(figsize=(10, 6))
    counts, bins, patches = ax.hist(
        healpix_gdf["cvh"], bins=50, alpha=0.7, edgecolor="black"
    )
    vmin, vcenter, vmax = 0.90, 1.00, 1.10
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i + 1]) / 2
        color = plt.cm.viridis(norm(bin_center))
        patch.set_facecolor(color)

    ax.axvline(x=1, color="red", linestyle="--", linewidth=2, label="Ideal (cvh = 1)")
    stats_text = (
        f"Mean: {healpix_gdf['cvh'].mean():.6f}\n"
        f"Std: {healpix_gdf['cvh'].std():.6f}\n"
        f"Min: {healpix_gdf['cvh'].min():.6f}\n"
        f"Max: {healpix_gdf['cvh'].max():.6f}"
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
    ax.set_xlabel("HEALPix CVH Compactness", fontsize=14)
    ax.set_ylabel("Number of cells", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()


def healpixinspect_cli():
    """Command-line interface for HEALPix cell inspection."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-r", "--resolution", dest="resolution", type=int, default=0)
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
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
        help="Antimeridian fixing method",
    )
    add_verbose_argument(parser)
    args = parser.parse_args()
    print(healpixinspect(args.resolution, fix_antimeridian=args.fix_antimeridian, verbose=args.verbose))


if __name__ == "__main__":
    healpixstats_cli()

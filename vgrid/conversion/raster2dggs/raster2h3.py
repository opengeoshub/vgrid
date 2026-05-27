"""
Raster to H3 Module

Convert raster data to H3 DGGS using either:

- **binning** — pixel centroids are aggregated into H3 cells (``stats`` required).
- **nearest_neighbour** — H3 grid over the raster bbox; each cell takes the value of
  the nearest raster pixel center (LAEA planar distance).
"""

from __future__ import annotations

import os
import argparse
from math import cos, radians

from tqdm import tqdm
import h3
import geopandas as gpd
import rasterio
from pyproj import datadir

from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    geodesic_dggs_to_geoseries,
    nearest_neighbour_from_grid,
)
from vgrid.utils.io import (
    validate_h3_resolution,
    convert_to_output_format,
    validate_raster_stats_option,
    normalize_raster2dggs_method,
    finalize_dggs_band_values,
)
from vgrid.utils.constants import (
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    DGGS_TYPES,
    MIN_CELL_AREA,
    RASTER_STATS_OPTIONS,
    RASTER2DGGS_METHODS,
)
from vgrid.conversion.dggs2geo.h32geo import h32geo
from vgrid.conversion.dggsresample.dggsresample import generate_grid

os.environ["PROJ_LIB"] = datadir.get_data_dir()

min_res = DGGS_TYPES["h3"]["min_res"]
max_res = DGGS_TYPES["h3"]["max_res"]


def get_nearest_h3_resolution(raster_path):
    """
    Automatically determine the optimal H3 resolution for a given raster.

    Analyzes the raster's pixel size and determines the most appropriate H3 resolution
    that best matches the raster's spatial resolution.
    """
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. H3 conversion requires a valid CRS."
            )
        pixel_width = transform.a
        pixel_height = -transform.e
        cell_size = pixel_width * pixel_height

        if crs.is_geographic:
            center_latitude = (src.bounds.top + src.bounds.bottom) / 2
            meter_per_degree_lat = 111_320
            meter_per_degree_lon = meter_per_degree_lat * cos(radians(center_latitude))
            pixel_width_m = pixel_width * meter_per_degree_lon
            pixel_height_m = pixel_height * meter_per_degree_lat
            cell_size = pixel_width_m * pixel_height_m

    min_diff = float("inf")
    nearest_resolution = min_res

    for res in range(min_res, max_res + 1):
        avg_area = h3.average_hexagon_area(res, unit="m^2")
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2h3_nearest_neighbour(
    raster_path: str,
    resolution: int,
    fix_antimeridian=None,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    h3_gdf = generate_grid(
        footprint, "h3", resolution, fix_antimeridian=fix_antimeridian
    )
    return nearest_neighbour_from_grid(raster_path, h3_gdf)


def _raster2h3_binning(
    raster_path: str,
    resolution: int,
    stats: str,
    fix_antimeridian=None,
) -> gpd.GeoDataFrame:
    """Bin pixel centroids into H3 cells and aggregate band values with ``stats``."""

    def cell_id(lat: float, lon: float):
        return h3.latlng_to_cell(lat, lon, resolution)

    h3_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to H3"
    )

    properties = []
    for h3_id, acc in tqdm(
        h3_acc.items(), desc="Converting raster to H3", unit=" cells"
    ):
        cell_polygon = h32geo(h3_id, fix_antimeridian=fix_antimeridian)
        num_edges = 6
        if h3.is_pentagon(h3_id):
            num_edges = 5
        base_props = geodesic_dggs_to_geoseries(
            "h3", h3_id, resolution, cell_polygon, num_edges
        )
        band_values = finalize_dggs_band_values(acc, stats)
        band_properties = {f"band_{i + 1}": band_values[i] for i in range(band_count)}
        base_props.update(band_properties)
        properties.append(base_props)

    if not properties:
        return gpd.GeoDataFrame(columns=["h3", "geometry"], crs="EPSG:4326")

    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2h3(
    raster_path,
    resolution=None,
    output_format="gpd",
    fix_antimeridian=None,
    method="binning",
    stats="mean",
):
    """
    Convert raster data to H3 DGGS format.

    Parameters
    ----------
    raster_path : str
        Path to the raster file (must have a defined CRS).
    resolution : int, optional
        H3 resolution [0..15]. If None, matched to pixel size.
    output_format : str, optional
        See :func:`~vgrid.utils.io.convert_to_output_format`.
    fix_antimeridian : str, optional
        Passed to H3 grid / geometry builders.
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"``: per-cell aggregation when multiple pixels map
        to one H3 cell (see ``RASTER_STATS_OPTIONS``). Ignored for ``nearest_neighbour``.
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_h3_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"H3 resolution: {resolution}")
    else:
        resolution = validate_h3_resolution(resolution)

    print(f"Method: {method}")

    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2h3_binning(
            raster_path, resolution, stats, fix_antimeridian=fix_antimeridian
        )
    else:
        gdf = _raster2h3_nearest_neighbour(
            raster_path, resolution, fix_antimeridian=fix_antimeridian
        )

    if gdf.empty:
        raise ValueError("No H3 cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2h3" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2h3_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to H3 DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"H3 resolution [{min_res}..{max_res}]",
    )

    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=RASTER2DGGS_METHODS,
        default="binning",
        help="binning: aggregate pixels into H3 cells; nearest_neighbour: nearest pixel center",
    )

    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
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
    parser.add_argument(
        "-s",
        "--stats",
        type=str,
        choices=RASTER_STATS_OPTIONS,
        default="mean",
        help="Band statistic for binning method only",
    )

    args = parser.parse_args()
    if not os.path.exists(args.raster):
        print(f"Error: The file {args.raster} does not exist.")
        return

    result = raster2h3(
        args.raster,
        args.resolution,
        args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        method=args.method,
        stats=args.stats,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    raster2h3_cli()

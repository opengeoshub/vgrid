"""
Raster to Maidenhead Module

Convert raster pixels to Maidenhead locator cells using either binning or
nearest-neighbour assignment.
"""

import argparse
import os
from math import cos, radians

import geopandas as gpd
import rasterio
from pyproj import datadir
from tqdm import tqdm

os.environ["PROJ_LIB"] = datadir.get_data_dir()

from vgrid.conversion.dggs2geo.maidenhead2geo import maidenhead2geo
from vgrid.conversion.latlon2dggs import latlon2maidenhead
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    nearest_neighbour_from_grid,
)
from vgrid.generator.maidenheadgrid import maidenhead_grid_within_bbox
from vgrid.stats.maidenheadstats import maidenhead_metrics
from vgrid.utils.constants import (
    DGGS_TYPES,
    MIN_CELL_AREA,
    OUTPUT_FORMATS,
    RASTER_STATS_OPTIONS,
    RASTER2DGGS_METHODS,
    STRUCTURED_FORMATS,
)
from vgrid.utils.io import (
    add_verbose_argument,
    convert_to_output_format,
    finalize_dggs_band_values,
    normalize_raster2dggs_method,
    validate_maidenhead_resolution,
    validate_raster_stats_option,
)

min_res = DGGS_TYPES["maidenhead"]["min_res"]
max_res = DGGS_TYPES["maidenhead"]["max_res"]


def get_nearest_maidenhead_resolution(raster_path):
    """Pick a Maidenhead resolution whose typical cell area is closest to the raster pixel area (m²)."""
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. Maidenhead conversion requires a valid CRS."
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
        _, _, avg_area, _ = maidenhead_metrics(res, unit="m")
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res
    return cell_size, nearest_resolution


def _raster2maidenhead_nearest_neighbour(
    raster_path: str, resolution: int,
    verbose=True,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    bbox = list(footprint.total_bounds)
    grid_gdf = maidenhead_grid_within_bbox(resolution, bbox, verbose=verbose)
    return nearest_neighbour_from_grid(raster_path, grid_gdf, verbose=verbose)


def _raster2maidenhead_binning(
    raster_path: str, resolution: int, stats: str,
    verbose=True,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        return latlon2maidenhead(lat, lon, resolution)

    maidenhead_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to Maidenhead", verbose=verbose
    )

    properties = []
    for maidenhead_id, acc in tqdm(
        maidenhead_acc.items(),
        desc="Converting raster to Maidenhead",
        unit=" cells",
        disable=not verbose,
    ):
        cell_polygon = maidenhead2geo(maidenhead_id)
        base_props = {"maidenhead": maidenhead_id, "geometry": cell_polygon}
        band_values = finalize_dggs_band_values(acc, stats)
        band_props = {f"band_{i + 1}": band_values[i] for i in range(band_count)}
        base_props.update(band_props)
        properties.append(base_props)

    if not properties:
        return gpd.GeoDataFrame(columns=["maidenhead", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2maidenhead(
    raster_path,
    resolution=None,
    output_format="gpd",
    method="binning",
    stats="mean",
    verbose=True,
):
    """
    Convert raster data to Maidenhead DGGS format.

    Parameters
    ----------
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"`` (see ``RASTER_STATS_OPTIONS``).
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_maidenhead_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest Maidenhead resolution determined: {resolution}")
    else:
        resolution = validate_maidenhead_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2maidenhead_binning(raster_path, resolution, stats, verbose=verbose)
    else:
        gdf = _raster2maidenhead_nearest_neighbour(raster_path, resolution, verbose=verbose)

    if gdf.empty:
        raise ValueError("No Maidenhead cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2maidenhead" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2maidenhead_cli():
    parser = argparse.ArgumentParser(
        description="Convert raster to Maidenhead DGGS (EPSG:4326 cells)"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"Maidenhead resolution [{min_res}..{max_res}]. Omit to infer from pixel size.",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=RASTER2DGGS_METHODS,
        default="binning",
        help="binning: aggregate pixels into Maidenhead cells; nearest_neighbour: nearest pixel center",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    parser.add_argument(
        "-s",
        "--stats",
        type=str,
        choices=RASTER_STATS_OPTIONS,
        default="mean",
        help="Band statistic for binning method only",
    )

    add_verbose_argument(parser)
    args = parser.parse_args()
    if not os.path.exists(args.raster):
        print(f"Error: The file {args.raster} does not exist.")
        return

    result = raster2maidenhead(
        args.raster,
        args.resolution,
        args.output_format,
        method=args.method,
        stats=args.stats,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    raster2maidenhead_cli()

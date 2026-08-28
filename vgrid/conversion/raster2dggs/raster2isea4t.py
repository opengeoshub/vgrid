"""
Raster to ISEA4T Module

Convert raster data to ISEA4T DGGS using either binning or nearest-neighbour assignment.

Note: This module is only supported on Windows systems due to OpenEaggr dependency.
"""

import os
import argparse
import platform
from math import cos, radians

from tqdm import tqdm
import geopandas as gpd
import rasterio
from pyproj import datadir

from vgrid.stats.isea4tstats import isea4t_metrics
from vgrid.utils.constants import ISEA4T_RES_ACCURACY_DICT
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    geodesic_dggs_metrics,
    nearest_neighbour_from_grid,
)
from vgrid.conversion.dggs2geo.isea4t2geo import isea4t2geo
from vgrid.conversion.dggsresample.dggsresample import generate_grid
from vgrid.utils.io import (
    add_verbose_argument,
    validate_isea4t_resolution,
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

os.environ["PROJ_LIB"] = datadir.get_data_dir()

min_res = DGGS_TYPES["isea4t"]["min_res"]
max_res = DGGS_TYPES["isea4t"]["max_res"]

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.enums.model import Model
    from vgrid.dggs.eaggr.shapes.lat_long_point import LatLongPoint

    isea4t_dggs = Eaggr(Model.ISEA4T)


def _require_windows():
    if platform.system() != "Windows":
        raise ValueError("ISEA4T is only supported on Windows systems.")


def get_nearest_isea4t_resolution(raster_path):
    """Automatically determine the optimal ISEA4T resolution for a given raster."""
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. ISEA4T conversion requires a valid CRS."
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
        _, _, avg_area, _ = isea4t_metrics(res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _isea4t_cell_id(lat, lon, resolution):
    max_accuracy = ISEA4T_RES_ACCURACY_DICT[39]
    lat_long_point = LatLongPoint(lat, lon, max_accuracy)
    isea4t_cell_max_accuracy = isea4t_dggs.convert_point_to_dggs_cell(lat_long_point)
    cell_id_len = resolution + 2
    isea4t_cell = DggsCell(isea4t_cell_max_accuracy._cell_id[:cell_id_len])
    return isea4t_cell._cell_id


def _raster2isea4t_nearest_neighbour(
    raster_path: str,
    resolution: int,
    fix_antimeridian=None,
    verbose=True,
) -> gpd.GeoDataFrame:
    _require_windows()
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(
        footprint, "isea4t", resolution, fix_antimeridian=fix_antimeridian, verbose=verbose
    )
    return nearest_neighbour_from_grid(raster_path, grid_gdf, verbose=verbose)


def _raster2isea4t_binning(
    raster_path: str,
    resolution: int,
    stats: str,
    fix_antimeridian=None,
    verbose=True,
) -> gpd.GeoDataFrame:
    _require_windows()

    def cell_id(lat, lon):
        return _isea4t_cell_id(lat, lon, resolution)

    isea4t_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to ISEA4T", verbose=verbose
    )

    properties = []
    for isea4t_id, acc in tqdm(
        isea4t_acc.items(),
        desc="Converting raster to ISEA4T",
        unit=" cells",
        disable=not verbose,
    ):
        cell_polygon = isea4t2geo(isea4t_id, fix_antimeridian=fix_antimeridian)
        num_edges = 3
        centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
            geodesic_dggs_metrics(cell_polygon, num_edges)
        )
        base_props = {
            "isea4t": isea4t_id,
            "resolution": resolution,
            "center_lat": centroid_lat,
            "center_lon": centroid_lon,
            "avg_edge_len": avg_edge_len,
            "cell_area": cell_area,
            "cell_perimeter": cell_perimeter,
            "geometry": cell_polygon,
        }
        band_values = finalize_dggs_band_values(acc, stats)
        band_properties = {f"band_{i + 1}": band_values[i] for i in range(band_count)}
        base_props.update(band_properties)
        properties.append(base_props)

    if not properties:
        return gpd.GeoDataFrame(columns=["isea4t", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2isea4t(
    raster_path,
    resolution=None,
    output_format="gpd",
    fix_antimeridian=None,
    method="binning",
    stats="mean",
    verbose=True,
):
    """
    Convert raster data to ISEA4T DGGS format.

    Parameters
    ----------
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"`` (see ``RASTER_STATS_OPTIONS``).
    """
    _require_windows()
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_isea4t_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest ISEA4T resolution determined: {resolution}")
    else:
        resolution = validate_isea4t_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2isea4t_binning(
            raster_path, resolution, stats, fix_antimeridian=fix_antimeridian, verbose=verbose
        )
    else:
        gdf = _raster2isea4t_nearest_neighbour(
            raster_path, resolution, fix_antimeridian=fix_antimeridian, verbose=verbose
        )

    if gdf.empty:
        raise ValueError("No ISEA4T cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2isea4t" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2isea4t_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to OpenEaggr ISEA4T DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"ISEA4T resolution [{min_res}..{max_res}]",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=RASTER2DGGS_METHODS,
        default="binning",
        help="binning: aggregate pixels into ISEA4T cells; nearest_neighbour: nearest pixel center",
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

    add_verbose_argument(parser)
    args = parser.parse_args()
    if not os.path.exists(args.raster):
        print(f"Error: The file {args.raster} does not exist.")
        return

    if platform.system() != "Windows":
        print("ISEA4T is only supported on Windows systems")
        return

    result = raster2isea4t(
        args.raster,
        args.resolution,
        args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        method=args.method,
        stats=args.stats,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    raster2isea4t_cli()

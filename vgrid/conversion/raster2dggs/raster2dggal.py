"""
Raster to DGGAL Module

Convert raster data to DGGAL DGGS using either binning or nearest-neighbour assignment.
"""

import os
import argparse
import math

from tqdm import tqdm
import geopandas as gpd
import rasterio
from pyproj import datadir
from dggal import *

from vgrid.conversion.latlon2dggs import latlon2dggal
from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo
from vgrid.conversion.dggsresample.dggsresample import generate_grid
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    geodesic_dggs_metrics,
    nearest_neighbour_from_grid,
)
from vgrid.stats.dggalstats import dggal_metrics
from vgrid.utils.io import (
    add_verbose_argument,
    convert_to_output_format,
    validate_dggal_resolution,
    validate_dggal_type,
    validate_raster_stats_option,
    normalize_raster2dggs_method,
    finalize_dggs_band_values,
)
from vgrid.utils.constants import (
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    DGGAL_TYPES,
    MIN_CELL_AREA,
    RASTER_STATS_OPTIONS,
    RASTER2DGGS_METHODS,
)

os.environ["PROJ_LIB"] = datadir.get_data_dir()

app = Application(appGlobals=globals())
pydggal_setup(app)


def get_nearest_dggal_resolution(dggs_type, raster_path):
    """Automatically determine the optimal DGGAL resolution for a given raster."""
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. DGGAL conversion requires a valid CRS."
            )
        pixel_width = transform.a
        pixel_height = -transform.e
        cell_size = pixel_width * pixel_height

        if crs.is_geographic:
            center_latitude = (src.bounds.top + src.bounds.bottom) / 2
            meter_per_degree_lat = 111_320
            meter_per_degree_lon = meter_per_degree_lat * math.cos(
                math.radians(center_latitude)
            )
            pixel_width_m = pixel_width * meter_per_degree_lon
            pixel_height_m = pixel_height * meter_per_degree_lat
            cell_size = pixel_width_m * pixel_height_m

    min_diff = float("inf")
    min_res = int(DGGAL_TYPES[dggs_type]["min_res"])
    max_res = int(DGGAL_TYPES[dggs_type]["max_res"])
    nearest_resolution = min_res
    for res in range(min_res, max_res + 1):
        _, _, avg_area, _ = dggal_metrics(dggs_type, res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = math.fabs(avg_area - cell_size)
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2dggal_nearest_neighbour(
    dggs_type: str,
    raster_path: str,
    resolution: int,
    split_antimeridian: bool = False,
    verbose=True,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(
        footprint,
        f"dggal_{dggs_type}",
        resolution,
        split_antimeridian=split_antimeridian, verbose=verbose
    )
    return nearest_neighbour_from_grid(raster_path, grid_gdf, verbose=verbose)


def _raster2dggal_binning(
    dggs_type: str,
    raster_path: str,
    resolution: int,
    stats: str,
    split_antimeridian: bool = False,
    verbose=True,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        try:
            return latlon2dggal(dggs_type, lat, lon, resolution)
        except Exception:
            return None

    zone_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to DGGAL", verbose=verbose
    )

    properties = []
    for zone_id, acc in tqdm(
        zone_acc.items(), desc="Converting raster to DGGAL", unit=" cells",
        disable=not verbose,
    ):
        try:
            dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
            dggrs = globals()[dggs_class_name]()
            zone = dggrs.getZoneFromTextID(zone_id)
            num_edges = dggrs.countZoneEdges(zone)
            cell_polygon = dggal2geo(
                dggs_type, zone_id, split_antimeridian=split_antimeridian
            )
            centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
                geodesic_dggs_metrics(cell_polygon, num_edges)
            )
            base_props = {
                f"dggal_{dggs_type}": zone_id,
                "resolution": resolution,
                "center_lat": centroid_lat,
                "center_lon": centroid_lon,
                "avg_edge_len": avg_edge_len,
                "cell_area": cell_area,
                "cell_perimeter": cell_perimeter,
                "geometry": cell_polygon,
            }
            band_values = finalize_dggs_band_values(acc, stats)
            band_properties = {
                f"band_{i + 1}": band_values[i] for i in range(band_count)
            }
            base_props.update(band_properties)
            properties.append(base_props)
        except Exception:
            continue

    if not properties:
        col = f"dggal_{dggs_type}"
        return gpd.GeoDataFrame(columns=[col, "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2dggal(
    dggs_type: str,
    raster_path,
    resolution: int | None = None,
    output_format: str = "gpd",
    split_antimeridian: bool = False,
    method: str = "binning",
    stats: str = "mean",
    verbose=True,
):
    """
    Convert raster data to DGGAL DGGS format.

    Parameters
    ----------
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"`` (see ``RASTER_STATS_OPTIONS``).
    """
    dggs_type = validate_dggal_type(dggs_type)
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_dggal_resolution(dggs_type, raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest {dggs_type.upper()} resolution determined: {resolution}")
    else:
        resolution = validate_dggal_resolution(dggs_type, resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2dggal_binning(
            dggs_type,
            raster_path,
            resolution,
            stats,
            split_antimeridian=split_antimeridian, verbose=verbose
        )
    else:
        gdf = _raster2dggal_nearest_neighbour(
            dggs_type,
            raster_path,
            resolution,
            split_antimeridian=split_antimeridian, verbose=verbose
        )

    if gdf.empty:
        raise ValueError("No DGGAL cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2dggal" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2dggal_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to DGGAL DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL type",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help="Resolution (integer). If omitted, auto-selected",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=RASTER2DGGS_METHODS,
        default="binning",
        help="binning: aggregate pixels into DGGAL cells; nearest_neighbour: nearest pixel center",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
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

    result = raster2dggal(
        args.dggs_type,
        args.raster,
        args.resolution,
        args.output_format,
        split_antimeridian=args.split_antimeridian,
        method=args.method,
        stats=args.stats,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    raster2dggal_cli()

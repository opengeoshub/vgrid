"""
Raster to DGGRID Module

Convert raster data to DGGRID DGGS using either binning or nearest-neighbour assignment.
"""

import os
import argparse
import json
import math

from tqdm import tqdm
import geopandas as gpd
import rasterio
from pyproj import datadir

from vgrid.conversion.dggsresample.dggsresample import generate_grid
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    dggrid_num_edges,
    footprint_gdf_from_raster,
    geodesic_dggs_metrics,
    nearest_neighbour_from_grid,
)
from vgrid.utils.io import (
    convert_to_output_format,
    validate_dggrid_type,
    validate_dggrid_resolution,
    create_dggrid_instance,
    validate_raster_stats_option,
    normalize_raster2dggs_method,
    finalize_dggs_band_values,
)
from vgrid.utils.constants import (
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    DGGRID_TYPES,
    MIN_CELL_AREA,
    RASTER_STATS_OPTIONS,
    RASTER2DGGS_METHODS,
)
from vgrid.stats.dggridstats import dggridstats
from vgrid.conversion.latlon2dggs import latlon2dggrid
from vgrid.conversion.dggs2geo.dggrid2geo import dggrid2geo

os.environ["PROJ_LIB"] = datadir.get_data_dir()


def get_nearest_dggrid_resolution(dggrid_instance, dggs_type, raster_path):
    """Automatically determine the optimal DGGRID resolution for a given raster."""
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. DGGRID conversion requires a valid CRS."
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
    min_res = int(DGGRID_TYPES[dggs_type]["min_res"])
    max_res = int(DGGRID_TYPES[dggs_type]["max_res"])
    nearest_resolution = min_res

    try:
        grid_stats = dggridstats(dggrid_instance, dggs_type, unit="m")
        for res in range(min_res, max_res + 1):
            res_stats = grid_stats[grid_stats["resolution"] == res]
            if res_stats.empty:
                continue
            avg_area_m2 = res_stats["area_m2"].iloc[0]
            if avg_area_m2 < MIN_CELL_AREA:
                break
            diff = math.fabs(avg_area_m2 - cell_size)
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res
    except Exception:
        nearest_resolution = min_res

    return cell_size, nearest_resolution


def _raster2dggrid_nearest_neighbour(
    dggrid_instance,
    dggs_type: str,
    raster_path: str,
    resolution: int,
    split_antimeridian: bool = False,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(
        footprint,
        f"dggrid_{dggs_type}",
        resolution,
        split_antimeridian=split_antimeridian,
    )
    return nearest_neighbour_from_grid(raster_path, grid_gdf)


def _raster2dggrid_binning(
    dggrid_instance,
    dggs_type: str,
    raster_path: str,
    resolution: int,
    stats: str,
    split_antimeridian: bool = False,
    aggregate: bool = False,
    options=None,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        try:
            return latlon2dggrid(
                dggrid_instance, dggs_type, lat, lon, resolution
            )
        except Exception:
            return None

    dggrid_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to DGGRID"
    )

    properties = []
    for dggrid_id, acc in tqdm(
        dggrid_acc.items(),
        desc="Converting raster to DGGRID",
        unit=" cells",
    ):
        try:
            cell_polygon = dggrid2geo(
                dggrid_instance,
                dggs_type,
                dggrid_id,
                resolution,
                split_antimeridian=split_antimeridian,
                aggregate=aggregate,
                options=options,
            )
            if not isinstance(cell_polygon, gpd.GeoDataFrame) or cell_polygon.empty:
                continue
            cell_geom = cell_polygon.iloc[0].geometry
            centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
                geodesic_dggs_metrics(cell_geom, dggrid_num_edges(dggs_type))
            )
            base_props = {
                f"dggrid_{dggs_type}": dggrid_id,
                "resolution": resolution,
                "center_lat": centroid_lat,
                "center_lon": centroid_lon,
                "avg_edge_len": avg_edge_len,
                "cell_area": cell_area,
                "cell_perimeter": cell_perimeter,
                "geometry": cell_geom,
            }
            band_values = finalize_dggs_band_values(acc, stats)
            band_properties = {
                f"band_{i + 1}": band_values[i] for i in range(band_count)
            }
            base_props.update(band_properties)
            properties.append(base_props)
        except Exception:
            continue

    col = f"dggrid_{dggs_type}"
    if not properties:
        return gpd.GeoDataFrame(columns=[col, "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2dggrid(
    dggrid_instance,
    dggs_type: str,
    raster_path,
    resolution: int | None = None,
    output_format: str = "gpd",
    split_antimeridian: bool = False,
    aggregate: bool = False,
    options=None,
    method: str = "binning",
    stats: str = "mean",
):
    """
    Convert raster data to DGGRID DGGS format.

    Parameters
    ----------
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"`` (see ``RASTER_STATS_OPTIONS``).
    """
    dggs_type = validate_dggrid_type(dggs_type)
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_dggrid_resolution(
            dggrid_instance, dggs_type, raster_path
        )
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest {dggs_type.upper()} resolution determined: {resolution}")
    else:
        resolution = validate_dggrid_resolution(dggs_type, resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2dggrid_binning(
            dggrid_instance,
            dggs_type,
            raster_path,
            resolution,
            stats,
            split_antimeridian=split_antimeridian,
            aggregate=aggregate,
            options=options,
        )
    else:
        gdf = _raster2dggrid_nearest_neighbour(
            dggrid_instance,
            dggs_type,
            raster_path,
            resolution,
            split_antimeridian=split_antimeridian,
        )

    if gdf.empty:
        raise ValueError("No DGGRID cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2dggrid" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2dggrid_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to DGGRID DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGRID_TYPES.keys(),
        help="DGGRID type",
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
        help="binning: aggregate pixels into DGGRID cells; nearest_neighbour: nearest pixel center",
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
        "-aggregate",
        "--aggregate",
        action="store_true",
        help="Aggregate the resulting polygons",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_from_cellids. "
        'Example: \'{"densification": 2}\'',
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

    dggrid_instance = create_dggrid_instance()

    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    result = raster2dggrid(
        dggrid_instance,
        args.dggs_type,
        args.raster,
        args.resolution,
        args.output_format,
        args.split_antimeridian,
        args.aggregate,
        options,
        method=args.method,
        stats=args.stats,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    raster2dggrid_cli()

"""
Raster to HEALPix Module

Convert raster data to HEALPix (UNIQ) DGGS using either:

- **binning** — pixel centroids are aggregated into HEALPix cells (``stats`` required).
- **nearest_neighbour** — HEALPix grid over the raster bbox; each cell takes the value
  of the nearest raster pixel center.
"""

from __future__ import annotations

import os
import argparse
from math import cos, radians

from tqdm import tqdm
import geopandas as gpd
import rasterio
from pyproj import datadir

from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    geodesic_dggs_metrics,
    nearest_neighbour_from_grid,
)
from vgrid.utils.io import (
    validate_healpix_resolution,
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
from vgrid.conversion.latlon2dggs import latlon2healpix
from vgrid.conversion.dggs2geo.healpix2geo import healpix2geo
from vgrid.conversion.dggsresample.dggsresample import generate_grid
from vgrid.dggs.healpix import nside2pixarea, order2nside, uniq2orderpix

os.environ["PROJ_LIB"] = datadir.get_data_dir()

min_res = DGGS_TYPES["healpix"]["min_res"]
max_res = DGGS_TYPES["healpix"]["max_res"]
EARTH_RADIUS_M = 6371008.8


def _healpix_avg_area_m2(resolution: int) -> float:
    """Average HEALPix cell area in square meters at a given order."""
    return nside2pixarea(order2nside(resolution)) * (EARTH_RADIUS_M**2)


def get_nearest_healpix_resolution(raster_path):
    """
    Automatically determine the optimal HEALPix resolution for a given raster.

    Returns
    -------
    tuple
        ``(cell_size_m2, resolution)``
    """
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. HEALPix conversion requires a valid CRS."
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
        avg_area = _healpix_avg_area_m2(res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2healpix_nearest_neighbour(
    raster_path: str,
    resolution: int,
    fix_antimeridian=None,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(
        footprint, "healpix", resolution, fix_antimeridian=fix_antimeridian
    )
    return nearest_neighbour_from_grid(raster_path, grid_gdf)


def _raster2healpix_binning(
    raster_path: str,
    resolution: int,
    stats: str,
    fix_antimeridian=None,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        return latlon2healpix(lat, lon, resolution)

    healpix_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to HEALPix"
    )

    properties = []
    for healpix_id, acc in tqdm(
        healpix_acc.items(),
        desc="Converting raster to HEALPix",
        unit=" cells",
    ):
        try:
            cell_polygon = healpix2geo(
                healpix_id, fix_antimeridian=fix_antimeridian
            )
            order, _ipix = uniq2orderpix(int(healpix_id))
            centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
                geodesic_dggs_metrics(cell_polygon, 4)
            )
            base_props = {
                "healpix": healpix_id,
                "resolution": order,
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
        return gpd.GeoDataFrame(columns=["healpix", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2healpix(
    raster_path,
    resolution=None,
    output_format="gpd",
    fix_antimeridian=None,
    method="binning",
    stats="mean",
):
    """
    Convert raster data to HEALPix UNIQ DGGS format.

    Parameters
    ----------
    raster_path : str
        Path to the raster file.
    resolution : int, optional
        HEALPix order [0..29]. If None, chosen to match raster pixel area.
    output_format : str, default "gpd"
        Output format.
    fix_antimeridian : str, optional
        Antimeridian fixing method.
    method : str, optional
        ``"binning"`` or ``"nearest_neighbour"``.
    stats : str, optional
        Aggregation for binning (see ``RASTER_STATS_OPTIONS``).
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_healpix_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest HEALPix resolution determined: {resolution}")
    else:
        resolution = validate_healpix_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2healpix_binning(
            raster_path, resolution, stats, fix_antimeridian=fix_antimeridian
        )
    else:
        gdf = _raster2healpix_nearest_neighbour(
            raster_path, resolution, fix_antimeridian=fix_antimeridian
        )

    if gdf.empty:
        raise ValueError("No HEALPix cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2healpix" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2healpix_cli():
    """Command-line interface for raster2healpix."""
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to HEALPix DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"HEALPix resolution [{min_res}..{max_res}]",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=RASTER2DGGS_METHODS,
        default="binning",
        help="binning: aggregate pixels into HEALPix cells; nearest_neighbour: nearest pixel center",
    )
    parser.add_argument(
        "-s",
        "--stats",
        type=str,
        choices=RASTER_STATS_OPTIONS,
        default="mean",
        help="Band statistic for binning method only",
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
    args = parser.parse_args()
    if not os.path.exists(args.raster):
        print(f"Error: The file {args.raster} does not exist.")
        return

    result = raster2healpix(
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
    raster2healpix_cli()

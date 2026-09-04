"""
Raster to DIGIPIN Module

Convert raster data to DIGIPIN DGGS using either:

- **binning** — pixel centroids aggregated into DIGIPIN cells (``stats`` required).
- **nearest_neighbour** — DIGIPIN grid over the raster bbox; each cell takes the
  nearest raster pixel center.
"""

import os
import argparse
from math import cos, radians

from tqdm import tqdm
import geopandas as gpd
import rasterio
from pyproj import datadir

from vgrid.conversion.latlon2dggs import latlon2digipin
from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    nearest_neighbour_from_grid,
)
from vgrid.generator.digipingrid import digipin_grid
from vgrid.stats.digipinstats import digipin_metrics
from vgrid.utils.io import (
    add_verbose_argument,
    validate_digipin_resolution,
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

min_res = DGGS_TYPES["digipin"]["min_res"]
max_res = DGGS_TYPES["digipin"]["max_res"]


def get_nearest_digipin_resolution(raster_path):
    """
    Automatically determine the optimal DIGIPIN resolution for a given raster.
    """
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. DIGIPIN conversion requires a valid CRS."
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
        _, _, avg_area, _ = digipin_metrics(res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2digipin_nearest_neighbour(
    raster_path: str, resolution: int,
    verbose=True,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    bbox = list(footprint.total_bounds)
    grid_gdf = digipin_grid(resolution, bbox, verbose=verbose)
    return nearest_neighbour_from_grid(raster_path, grid_gdf, verbose=verbose)


def _raster2digipin_binning(
    raster_path: str, resolution: int, stats: str,
    verbose=True,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        return latlon2digipin(lat, lon, resolution)

    digipin_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to DIGIPIN", verbose=verbose
    )

    properties = []
    for digipin_id, acc in tqdm(
        digipin_acc.items(),
        desc="Converting raster to DIGIPIN",
        unit=" cells",
        disable=not verbose,
    ):
        cell_polygon = digipin2geo(digipin_id)
        if isinstance(cell_polygon, str):
            continue
        base_props = {"digipin": digipin_id, "geometry": cell_polygon}
        band_values = finalize_dggs_band_values(acc, stats)
        band_props = {f"band_{i + 1}": band_values[i] for i in range(band_count)}
        base_props.update(band_props)
        properties.append(base_props)

    if not properties:
        return gpd.GeoDataFrame(columns=["digipin", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2digipin(
    raster_path,
    resolution=None,
    output_format="gpd",
    method="binning",
    stats="mean",
    verbose=True,
):
    """
    Convert raster data to DIGIPIN DGGS format.

    Parameters
    ----------
    raster_path : str
        Path to the raster file.
    resolution : int, optional
        DIGIPIN resolution [1..10]. If None, matched to pixel size.
    output_format : str, optional
        See :func:`~vgrid.utils.io.convert_to_output_format`.
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"`` (see ``RASTER_STATS_OPTIONS``).
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_digipin_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest DIGIPIN resolution determined: {resolution}")
    else:
        resolution = validate_digipin_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2digipin_binning(raster_path, resolution, stats, verbose=verbose)
    else:
        gdf = _raster2digipin_nearest_neighbour(raster_path, resolution, verbose=verbose)

    if gdf.empty:
        raise ValueError("No DIGIPIN cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2digipin" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2digipin_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to DIGIPIN DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"DIGIPIN resolution [{min_res}..{max_res}]",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=RASTER2DGGS_METHODS,
        default="binning",
        help="binning: aggregate pixels into DIGIPIN cells; nearest_neighbour: nearest pixel center",
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

    result = raster2digipin(
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
    raster2digipin_cli()

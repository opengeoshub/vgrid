"""
Raster to Tilecode Module

This module provides functionality to convert raster data to Tilecode DGGS format with automatic resolution determination and multi-band support.

Key Functions:
    raster2tilecode: Main conversion function with multiple output formats
    get_nearest_tilecode_resolution: Automatically determines optimal Tilecode resolution
    raster2tilecode_cli: Command-line interface for conversion process
"""

import os
import argparse
import re
from math import cos, radians
from tqdm import tqdm
from shapely.geometry import Polygon
import rasterio
import geopandas as gpd
from pyproj import datadir
from vgrid.dggs import tilecode, mercantile
from vgrid.stats.tilecodestats import tilecode_metrics
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    nearest_neighbour_from_grid,
)
from vgrid.conversion.dggsresample.dggsresample import generate_grid
from vgrid.utils.io import (
    add_verbose_argument,
    validate_tilecode_resolution,
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

min_res = DGGS_TYPES["tilecode"]["min_res"]
max_res = DGGS_TYPES["tilecode"]["max_res"]


def get_nearest_tilecode_resolution(raster_path):
    """
    Automatically determine the optimal Tilecode resolution for a given raster.

    Analyzes the raster's pixel size and determines the most appropriate Tilecode resolution
    that best matches the raster's spatial resolution.

    Parameters
    ----------
    raster_path : str
        Path to the raster file to analyze.

    Returns
    -------
    tuple
        A tuple containing (cell_size, resolution) where:
        - cell_size: The calculated cell size in square meters
        - resolution: The optimal Tilecode resolution level

    Examples
    --------
    >>> cell_size, resolution = get_nearest_tilecode_resolution("data.tif")
    >>> print(f"Cell size: {cell_size} m², Resolution: {resolution}")
    Cell size: 1000000.0 m², Resolution: 5
    """
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. Tilecode conversion requires a valid CRS."
            )
        pixel_width = transform.a
        pixel_height = -transform.e
        cell_size = pixel_width * pixel_height

        if crs.is_geographic:
            # Latitude of the raster center
            center_latitude = (src.bounds.top + src.bounds.bottom) / 2
            # Convert degrees to meters
            meter_per_degree_lat = 111_320  # Roughly 1 degree latitude in meters
            meter_per_degree_lon = meter_per_degree_lat * cos(radians(center_latitude))

            pixel_width_m = pixel_width * meter_per_degree_lon
            pixel_height_m = pixel_height * meter_per_degree_lat
            cell_size = pixel_width_m * pixel_height_m

    min_diff = float("inf")
    nearest_resolution = min_res

    for res in range(min_res, max_res + 1):
        _, _, avg_area, _ = tilecode_metrics(res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        # If the difference is smaller than the current minimum, update the nearest resolution
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2tilecode_nearest_neighbour(raster_path, resolution,
    verbose=True,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(footprint, "tilecode", resolution, verbose=verbose)
    return nearest_neighbour_from_grid(raster_path, grid_gdf, verbose=verbose)


def _raster2tilecode_binning(raster_path, resolution, stats,
    verbose=True,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        return tilecode.latlon2tilecode(lat, lon, resolution)

    tilecode_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to Tilecode", verbose=verbose
    )

    properties = []
    for tilecode_id, acc in tqdm(
        tilecode_acc.items(),
        desc="Converting raster to Tilecode",
        unit=" cells",
        disable=not verbose,
    ):
        match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
        if match:
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                cell_polygon = Polygon(
                    [
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat],
                    ]
                )
                base_props = {"tilecode": tilecode_id, "geometry": cell_polygon}
                band_values = finalize_dggs_band_values(acc, stats)
                band_props = {
                    f"band_{i + 1}": band_values[i] for i in range(band_count)
                }
                base_props.update(band_props)
                properties.append(base_props)

    if not properties:
        return gpd.GeoDataFrame(columns=["tilecode", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2tilecode(
    raster_path, resolution=None, output_format="gpd", method="binning", stats="mean",
    verbose=True,
):
    """
    Convert raster data to Tilecode DGGS format.

    Converts raster data to Tilecode DGGS format with automatic resolution
    determination and multi-band support. Each pixel is assigned to a Tilecode cell and
    the first sample value per cell is preserved.

    Parameters
    ----------
    raster_path : str
        Path to the raster file to convert.
    resolution : int, optional
        Tilecode resolution level. If None, automatically determined based on raster pixel size.
        Valid range: 0-26.
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"``: per-cell aggregation (see ``RASTER_STATS_OPTIONS``).

    Returns
    -------
    geopandas.GeoDataFrame or str or dict
        The converted data in the specified format. Each row represents a Tilecode cell
        with geometry and band values from the original raster.

    Examples
    --------
    >>> # Convert with automatic resolution
    >>> result = raster2tilecode("data.tif")
    >>> print(f"Converted {len(result)} Tilecode cells")

    >>> # Convert with specific resolution
    >>> result = raster2tilecode("data.tif", resolution=10)

    >>> # Convert to GeoJSON file
    >>> result = raster2tilecode("data.tif", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_tilecode_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest Tilecode resolution determined: {resolution}")
    else:
        resolution = validate_tilecode_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2tilecode_binning(raster_path, resolution, stats, verbose=verbose)
    else:
        gdf = _raster2tilecode_nearest_neighbour(raster_path, resolution, verbose=verbose)

    if gdf.empty:
        raise ValueError("No Tilecode cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2tilecode" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2tilecode_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to Tilecode DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"Tilecode resolution [{min_res}..{max_res}]",
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
        help="binning: aggregate pixels into Tilecode cells; nearest_neighbour: nearest pixel center",
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

    result = raster2tilecode(
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
    raster2tilecode_cli()

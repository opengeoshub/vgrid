"""
Raster to S2 Module

This module provides functionality to convert raster data to S2 DGGS format with automatic resolution determination and multi-band support.

Key Functions:
    raster2s2: Main conversion function with multiple output formats
    get_nearest_s2_resolution: Automatically determines optimal S2 resolution
    raster2s2_cli: Command-line interface for conversion process
"""

import os
import argparse
from tqdm import tqdm
import rasterio
from vgrid.dggs import s2
from vgrid.stats.s2stats import s2_metrics
from math import cos, radians
from vgrid.utils.io import (
    validate_s2_resolution,
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
import geopandas as gpd
from pyproj import datadir
from vgrid.conversion.dggs2geo.s22geo import s22geo
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    geodesic_dggs_metrics,
    nearest_neighbour_from_grid,
)
from vgrid.conversion.dggsresample.dggsresample import generate_grid

os.environ["PROJ_LIB"] = datadir.get_data_dir()
min_res = DGGS_TYPES["s2"]["min_res"]
max_res = DGGS_TYPES["s2"]["max_res"]


def get_nearest_s2_resolution(raster_path):
    """
    Automatically determine the optimal S2 resolution for a given raster.

    Analyzes the raster's pixel size and determines the most appropriate S2 resolution
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
        - resolution: The optimal S2 resolution level

    Examples
    --------
    >>> cell_size, resolution = get_nearest_s2_resolution("data.tif")
    >>> print(f"Cell size: {cell_size} m², Resolution: {resolution}")
    Cell size: 1000000.0 m², Resolution: 5
    """
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. S2 conversion requires a valid CRS."
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
    # Check resolutions from 0 to 15
    nearest_resolution = min_res

    for res in range(min_res, max_res + 1):
        _, _, avg_area, _ = s2_metrics(res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        # If the difference is smaller than the current minimum, update the nearest resolution
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2s2_nearest_neighbour(
    raster_path, resolution, fix_antimeridian=None
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(
        footprint, "s2", resolution, fix_antimeridian=fix_antimeridian
    )
    return nearest_neighbour_from_grid(raster_path, grid_gdf)


def _raster2s2_binning(
    raster_path, resolution, stats, fix_antimeridian=None
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        lat_lng = s2.LatLng.from_degrees(lat, lon)
        s2_id = s2.CellId.from_lat_lng(lat_lng).parent(resolution)
        return s2.CellId.to_token(s2_id)

    s2_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to S2"
    )

    properties = []
    for s2_token, acc in tqdm(
        s2_acc.items(), desc="Converting raster to S2", unit=" cells"
    ):
        cell_polygon = s22geo(s2_token, fix_antimeridian=fix_antimeridian)
        num_edges = 4
        centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
            geodesic_dggs_metrics(cell_polygon, num_edges)
        )
        base_props = {
            "s2": s2_token,
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
        return gpd.GeoDataFrame(columns=["s2", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2s2(
    raster_path,
    resolution=None,
    output_format="gpd",
    fix_antimeridian=None,
    method="binning",
    stats="mean",
):
    """
    Convert raster data to S2 DGGS format.

    Converts raster data to S2 DGGS format with automatic resolution
    determination and multi-band support. Each pixel is assigned to an S2 cell and
    the first sample value per cell is preserved.

    Parameters
    ----------
    raster_path : str
        Path to the raster file to convert.
    resolution : int, optional
        S2 resolution level. If None, automatically determined based on raster pixel size.
        Valid range: 0-30.
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"``: per-cell aggregation (see ``RASTER_STATS_OPTIONS``).
    Returns
    -------
    geopandas.GeoDataFrame or str or dict
        The converted data in the specified format. Each row represents an S2 cell
        with geometry and band values from the original raster.

    Examples
    --------
    >>> # Convert with automatic resolution
    >>> result = raster2s2("data.tif")
    >>> print(f"Converted {len(result)} S2 cells")

    >>> # Convert with specific resolution
    >>> result = raster2s2("data.tif", resolution=10)

    >>> # Convert to GeoJSON file
    >>> result = raster2s2("data.tif", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_s2_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest S2 resolution determined: {resolution}")
    else:
        resolution = validate_s2_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2s2_binning(
            raster_path, resolution, stats, fix_antimeridian=fix_antimeridian
        )
    else:
        gdf = _raster2s2_nearest_neighbour(
            raster_path, resolution, fix_antimeridian=fix_antimeridian
        )

    if gdf.empty:
        raise ValueError("No S2 cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2s2" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2s2_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to S2 DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"S2 resolution [{min_res}..{max_res}]",
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
        help="binning: aggregate pixels into S2 cells; nearest_neighbour: nearest pixel center",
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
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
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

    result = raster2s2(
        args.raster,
        args.resolution,
        args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        method=args.method,
        stats=args.stats,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    raster2s2_cli()

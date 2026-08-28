"""
Raster to RHEALPix Module

This module provides functionality to convert raster data to RHEALPix (Rectified HEALPix) DGGS format with automatic resolution determination and multi-band support.

Key Functions:
    raster2rhealpix: Main conversion function with multiple output formats
    get_nearest_rhealpix_resolution: Automatically determines optimal RHEALPix resolution
    raster2rhealpix_cli: Command-line interface for conversion process
"""

import os
import argparse
from tqdm import tqdm
from vgrid.stats.rhealpixstats import rhealpix_metrics
from vgrid.utils.io import (
    add_verbose_argument,
    validate_rhealpix_resolution,
    convert_to_output_format,
    validate_raster_stats_option,
    normalize_raster2dggs_method,
    finalize_dggs_band_values,
)
from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo
from vgrid.utils.constants import (
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    DGGS_TYPES,
    MIN_CELL_AREA,
    RASTER_STATS_OPTIONS,
    RASTER2DGGS_METHODS,
)
from math import cos, radians
from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.conversion.latlon2dggs import latlon2rhealpix
import geopandas as gpd
from pyproj import datadir
from vgrid.utils.geometry import (
    accumulate_raster_pixels,
    footprint_gdf_from_raster,
    geodesic_dggs_metrics,
    nearest_neighbour_from_grid,
)
from vgrid.conversion.dggsresample.dggsresample import generate_grid

os.environ["PROJ_LIB"] = datadir.get_data_dir()
import rasterio

E = WGS84_ELLIPSOID
rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
min_res = DGGS_TYPES["rhealpix"]["min_res"]
max_res = DGGS_TYPES["rhealpix"]["max_res"]


def get_nearest_rhealpix_resolution(raster_path):
    """
    Automatically determine the optimal RHEALPix resolution for a given raster.

    Analyzes the raster's pixel size and determines the most appropriate RHEALPix resolution
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
        - resolution: The optimal RHEALPix resolution level

    Examples
    --------
    >>> cell_size, resolution = get_nearest_rhealpix_resolution("data.tif")
    >>> print(f"Cell size: {cell_size} m², Resolution: {resolution}")
    Cell size: 1000000.0 m², Resolution: 5
    """
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        if crs is None:
            raise ValueError(
                "Raster CRS is undefined. rHEALPix conversion requires a valid CRS."
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
        _, _, avg_area, _ = rhealpix_metrics(res)
        if avg_area < MIN_CELL_AREA:
            break
        diff = abs(avg_area - cell_size)
        # If the difference is smaller than the current minimum, update the nearest resolution
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return cell_size, nearest_resolution


def _raster2rhealpix_nearest_neighbour(
    raster_path: str,
    resolution: int,
    fix_antimeridian=None,
    verbose=True,
) -> gpd.GeoDataFrame:
    footprint = footprint_gdf_from_raster(raster_path)
    grid_gdf = generate_grid(
        footprint, "rhealpix", resolution, fix_antimeridian=fix_antimeridian, verbose=verbose
    )
    return nearest_neighbour_from_grid(raster_path, grid_gdf, verbose=verbose)


def _raster2rhealpix_binning(
    raster_path: str,
    resolution: int,
    stats: str,
    fix_antimeridian=None,
    verbose=True,
) -> gpd.GeoDataFrame:
    def cell_id(lat, lon):
        return latlon2rhealpix(lat, lon, resolution)

    rhealpix_acc, band_count = accumulate_raster_pixels(
        raster_path, cell_id, stats, desc="Binning raster blocks to rHEALPix", verbose=verbose
    )

    properties = []
    for rhealpix_id, acc in tqdm(
        rhealpix_acc.items(),
        desc="Converting raster to rHEALPix",
        unit=" cells",
        disable=not verbose,
    ):
        cell_polygon = rhealpix2geo(rhealpix_id, fix_antimeridian=fix_antimeridian)
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        num_edges = 4
        if rhealpix_cell.ellipsoidal_shape() == "dart":
            num_edges = 3
        centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
            geodesic_dggs_metrics(cell_polygon, num_edges)
        )
        base_props = {
            "rhealpix": rhealpix_id,
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
        return gpd.GeoDataFrame(columns=["rhealpix", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(properties, geometry="geometry", crs="EPSG:4326")


def raster2rhealpix(
    raster_path,
    resolution=None,
    output_format="gpd",
    fix_antimeridian=None,
    method="binning",
    stats="mean",
    verbose=True,
):
    """
    Convert raster data to RHEALPix DGGS format.

    Converts raster data to RHEALPix (Rectified HEALPix) DGGS format with automatic resolution
    determination and multi-band support.

    Parameters
    ----------
    raster_path : str
        Path to the raster file to convert.
    resolution : int, optional
        RHEALPix resolution level. If None, automatically determined based on raster pixel size.
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
    method : str, optional
        ``"binning"`` (default) or ``"nearest_neighbour"`` (see ``RASTER2DGGS_METHODS``).
    stats : str, optional
        Used when ``method="binning"``: per-cell aggregation (see ``RASTER_STATS_OPTIONS``).

    Returns
    -------
    geopandas.GeoDataFrame or str or dict
        The converted data in the specified format. Each row represents a RHEALPix cell
        with geometry and band values from the original raster.

    Examples
    --------
    >>> # Convert with automatic resolution
    >>> result = raster2rhealpix("data.tif")
    >>> print(f"Converted {len(result)} RHEALPix cells")

    >>> # Convert with specific resolution
    >>> result = raster2rhealpix("data.tif", resolution=10)

    >>> # Convert to GeoJSON file
    >>> result = raster2rhealpix("data.tif", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    method = normalize_raster2dggs_method(method)

    if resolution is None:
        cell_size, resolution = get_nearest_rhealpix_resolution(raster_path)
        print(f"Cell size: {cell_size} m2")
        print(f"Nearest rHEALPix resolution determined: {resolution}")
    else:
        resolution = validate_rhealpix_resolution(resolution)

    print(f"Method: {method}")
    if method == "binning":
        stats = validate_raster_stats_option(stats)
        print(f"Stats: {stats}")
        gdf = _raster2rhealpix_binning(
            raster_path, resolution, stats, fix_antimeridian=fix_antimeridian, verbose=verbose
        )
    else:
        gdf = _raster2rhealpix_nearest_neighbour(
            raster_path, resolution, fix_antimeridian=fix_antimeridian, verbose=verbose
        )

    if gdf.empty:
        raise ValueError("No rHEALPix cells were produced from the raster.")

    base_name = os.path.splitext(os.path.basename(raster_path))[0]
    output_name = f"{base_name}2rhealpix" if output_format is not None else None
    return convert_to_output_format(gdf, output_format, output_name)


def raster2rhealpix_cli():
    """Command line interface for raster2rhealpix"""
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to rHEALPix DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help=f"rHEALPix resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="binning: aggregate pixels into rHEALPix cells; nearest_neighbour: nearest pixel center",
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
    add_verbose_argument(parser)
    args = parser.parse_args()
    if not os.path.exists(args.raster):
        raise FileNotFoundError(f"The file {args.raster} does not exist.")

    result = raster2rhealpix(
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
    raster2rhealpix_cli()

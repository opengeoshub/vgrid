"""
RHEALPix to Geometry Module

This module provides functionality to convert RHEALPix (Rectified HEALPix) cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    rhealpix2geo: Convert RHEALPix cell IDs to Shapely Polygons
    rhealpix2geojson: Convert RHEALPix cell IDs to GeoJSON FeatureCollection
    rhealpix2geo_cli: Command-line interface for polygon conversion
    rhealpix2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from vgrid.utils.geometry import dggs_geojson_feature, rhealpix_cell_to_polygon
from vgrid.utils.geometry import shift_balanced, shift_west, shift_east
from vgrid.utils.antimeridian import fix_polygon
from vgrid.utils.io import (
    add_rhealpix_n_side_argument,
    get_rhealpix_dggs,
    rhealpix_cell_from_id,
)
from pyproj import Geod
from vgrid.utils.constants import FIX_ANTIMERIDIAN_CHOICES

geod = Geod(ellps="WGS84")


def rhealpix2geo(rhealpix_ids, fix_antimeridian=None, N_side=3):
    """
    Convert RHEALPix cell IDs to Shapely geometry objects.

    Accepts a single rhealpix_id (string) or a list of rhealpix_ids. For each valid RHEALPix cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    rhealpix_ids : str or list of str
        RHEALPix cell ID(s) to convert. Can be a single string or a list of strings.
        Each ID should be a string starting with 'R' followed by numeric digits.
        Example format: "R31260335553825"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    N_side : int, default 3
        Children per cell edge (2 or 3). Must match the DGGS used to create the IDs.

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single RHEALPix cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding RHEALPix cell.

    Examples
    --------
    >>> rhealpix2geo("R31260335553825")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> rhealpix2geo(["R31260335553825", "R31260335553826"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(rhealpix_ids, str):
        rhealpix_ids = [rhealpix_ids]
    rhealpix_dggs = get_rhealpix_dggs(N_side=N_side)
    rhealpix_polygons = []
    for rhealpix_id in rhealpix_ids:
        try:
            rhealpix_cell = rhealpix_cell_from_id(rhealpix_id, dggs=rhealpix_dggs)
            cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
            if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
                cell_polygon = shift_balanced(
                    cell_polygon, threshold_west=-149, threshold_east=149
                )
            elif fix_antimeridian == "shift_west":
                cell_polygon = shift_west(cell_polygon, threshold=-149)
            elif fix_antimeridian == "shift_east":
                cell_polygon = shift_east(cell_polygon, threshold=149)
            elif fix_antimeridian == "split":
                cell_polygon = fix_polygon(cell_polygon)
            rhealpix_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(rhealpix_polygons) == 1:
        return rhealpix_polygons[0]
    return rhealpix_polygons


def rhealpix2geo_cli():
    """
    Command-line interface for converting RHEALPix cell IDs to Shapely Polygons.

    This function provides a command-line interface that accepts multiple RHEALPix
    cell IDs as command-line arguments and returns the corresponding Shapely
    Polygon objects.

    Returns:
        list: A list of Shapely Polygon objects representing the converted cells.

    Usage:
        rhealpix2geo R31260335553825 R31260335553826

    Note:
        This function is designed to be called from the command line and requires
        RHEALPix cell IDs as command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Convert Rhealpix to Geometry")
    parser.add_argument(
        "rhealpix",
        nargs="+",
        help="Input Rhealpix (string or list of strings)",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    add_rhealpix_n_side_argument(parser)
    args = parser.parse_args()
    polys = rhealpix2geo(args.rhealpix, args.fix_antimeridian, N_side=args.N_side)
    return polys


def rhealpix2geojson(rhealpix_ids, fix_antimeridian=None, N_side=3, cell_metrics=False):
    """
    Convert RHEALPix cell IDs to GeoJSON FeatureCollection.

    Accepts a single rhealpix_id (string) or a list of rhealpix_ids. For each valid
    RHEALPix cell ID, creates a GeoJSON feature with polygon geometry and
    cell metadata. Skips invalid or error-prone cells.

    Parameters
    ----------
    rhealpix_ids : str or list of str
        RHEALPix cell ID(s) to convert. Can be a single string or a list of strings.
        Each ID should be a string starting with 'R' followed by numeric digits.
        Example format: "R31260335553825"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    N_side : int, default 3
        Children per cell edge (2 or 3). Must match the DGGS used to create the IDs.

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing Polygon features for each valid
        RHEALPix cell. Each feature includes properties with cell ID, resolution,
        and other metadata.
    """
    if isinstance(rhealpix_ids, str):
        rhealpix_ids = [rhealpix_ids]
    rhealpix_features = []
    rhealpix_dggs = get_rhealpix_dggs(N_side=N_side)
    for rhealpix_id in rhealpix_ids:
        try:
            rhealpix_cell = rhealpix_cell_from_id(rhealpix_id, dggs=rhealpix_dggs)
            resolution = rhealpix_cell.resolution
            cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
            if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
                cell_polygon = shift_balanced(
                    cell_polygon, threshold_west=-128, threshold_east=160
                )
            elif fix_antimeridian == "shift_west":
                cell_polygon = shift_west(cell_polygon, threshold=-128)
            elif fix_antimeridian == "shift_east":
                cell_polygon = shift_east(cell_polygon, threshold=160)
            elif fix_antimeridian == "split":
                cell_polygon = fix_polygon(cell_polygon)
            num_edges = 4
            if rhealpix_cell.ellipsoidal_shape == "dart":
                num_edges = 3
            feature = dggs_geojson_feature(
                "rhealpix",
                rhealpix_id,
                resolution,
                cell_polygon,
                cell_metrics,
                num_edges,
            )
            rhealpix_features.append(feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": rhealpix_features}


def rhealpix2geojson_cli():
    """
    Command-line interface for converting RHEALPix cell IDs to GeoJSON FeatureCollection.

    This function provides a command-line interface that accepts multiple RHEALPix
    cell IDs as command-line arguments and returns the corresponding GeoJSON
    FeatureCollection as a JSON string.

    Returns:
        None: Prints the GeoJSON FeatureCollection to stdout.

    Usage:
        rhealpix2geojson R31260335553825 R31260335553826

    Note:
        This function is designed to be called from the command line and requires
        RHEALPix cell IDs as command-line arguments. The output is printed to
        stdout as a formatted JSON string.
    """
    parser = argparse.ArgumentParser(description="Convert Rhealpix to GeoJSON")
    parser.add_argument(
        "rhealpix",
        nargs="+",
        help="Input Rhealpix",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    add_rhealpix_n_side_argument(parser)
    parser.add_argument(
        "-cell_metrics",
        "--cell_metrics",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Include geodesic or graticule cell metrics. Default is off.",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(
        rhealpix2geojson(
            args.rhealpix,
            args.fix_antimeridian,
            N_side=args.N_side,
            cell_metrics=args.cell_metrics,
        )
    )
    print(geojson_data)

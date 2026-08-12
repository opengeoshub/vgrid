"""
HEALPix to Geometry Module

This module provides functionality to convert HEALPix UNIQ cell IDs to Shapely
Polygons and GeoJSON FeatureCollection.

Key Functions:
    healpix2geo: Convert HEALPix cell IDs to Shapely Polygons
    healpix2geojson: Convert HEALPix cell IDs to GeoJSON FeatureCollection
    healpix2geo_cli: Command-line interface for polygon conversion
    healpix2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from shapely.geometry import Polygon

from vgrid.dggs.healpix import (
    cornersNestLonLat,
    nside2npix,
    order2nside,
    uniq2orderpix,
)
from vgrid.utils.geometry import (
    geodesic_dggs_to_feature,
    shift_balanced,
    shift_west,
    shift_east,
)
from vgrid.utils.antimeridian import fix_polygon
from vgrid.utils.io import validate_healpix_resolution


def _normalize_ids(healpix_ids):
    """Normalize a single ID or sequence into a list of ints."""
    if isinstance(healpix_ids, (int, str)):
        healpix_ids = [healpix_ids]
    return [int(healpix_id) for healpix_id in healpix_ids]


def _decode_uniq(healpix_id):
    """Decode a UNIQ ID into (order, nested ipix, nside)."""
    order, ipix = uniq2orderpix(healpix_id)
    order = validate_healpix_resolution(order)
    nside = order2nside(order)
    if ipix < 0 or ipix >= nside2npix(nside):
        raise ValueError(f"Invalid nested pixel index {ipix} for order {order}")
    return order, ipix, nside


def _healpix_uniq_to_polygon(healpix_id):
    """Build a Shapely Polygon from a HEALPix UNIQ cell ID."""
    _order, ipix, nside = _decode_uniq(healpix_id)
    corners = cornersNestLonLat(nside, ipix)
    return Polygon([(lon, lat) for lon, lat in corners])


def _apply_antimeridian(cell_polygon, fix_antimeridian):
    """Apply optional antimeridian fixing to a cell polygon."""
    if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
        return shift_balanced(
            cell_polygon, threshold_west=-130, threshold_east=146
        )
    if fix_antimeridian == "shift_west":
        return shift_west(cell_polygon, threshold=-130)
    if fix_antimeridian == "shift_east":
        return shift_east(cell_polygon, threshold=146)
    if fix_antimeridian == "split":
        return fix_polygon(cell_polygon)
    return cell_polygon


def healpix2geo(healpix_ids, fix_antimeridian=None):
    """
    Convert HEALPix UNIQ cell IDs to Shapely geometry objects.

    Accepts a single healpix_id (int/str) or a list of IDs. Each UNIQ ID encodes
    both resolution and nested pixel index, so no separate resolution argument is
    required. Skips invalid or error-prone cells.

    Parameters
    ----------
    healpix_ids : int, str, or list of int/str
        HEALPix UNIQ cell ID(s) to convert (as returned by latlon2healpix).
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east,
        split, or none. Defaults to None (no fixing).

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single HEALPix cell ID is provided, returns a single Shapely Polygon.
        If a list of IDs is provided, returns a list of Shapely Polygons.

    Examples
    --------
    >>> healpix2geo(9941583)
    <shapely.geometry.polygon.Polygon object at ...>

    >>> healpix2geo([9941583, 9941584])
    [<shapely.geometry.polygon.Polygon object at ...>, ...]
    """
    single = isinstance(healpix_ids, (int, str))
    healpix_ids = _normalize_ids(healpix_ids)

    healpix_polygons = []
    for healpix_id in healpix_ids:
        try:
            cell_polygon = _healpix_uniq_to_polygon(healpix_id)
            cell_polygon = _apply_antimeridian(cell_polygon, fix_antimeridian)
            healpix_polygons.append(cell_polygon)
        except Exception:
            continue

    if single and len(healpix_polygons) == 1:
        return healpix_polygons[0]
    return healpix_polygons


def healpix2geo_cli():
    """
    Command-line interface for healpix2geo supporting multiple HEALPix cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert HEALPix UNIQ cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "healpix",
        nargs="+",
        help="Input HEALPix UNIQ cell ID(s), e.g., healpix2geo 9941583 9941584",
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
    args = parser.parse_args()
    polys = healpix2geo(args.healpix, fix_antimeridian=args.fix_antimeridian)
    return polys


def healpix2geojson(healpix_ids, fix_antimeridian=None):
    """
    Convert HEALPix UNIQ cell IDs to GeoJSON FeatureCollection.

    Accepts a single healpix_id (int/str) or a list of IDs. Each UNIQ ID encodes
    both resolution and nested pixel index. Skips invalid or error-prone cells.

    Parameters
    ----------
    healpix_ids : int, str, or list of int/str
        HEALPix UNIQ cell ID(s) to convert (as returned by latlon2healpix).
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east,
        split, or none. Defaults to None (no fixing).

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid cell.
        Each feature includes geometry and properties (cell ID, resolution, metrics).

    Examples
    --------
    >>> healpix2geojson(9941583)
    {'type': 'FeatureCollection', 'features': [...]}

    >>> healpix2geojson([9941583, 9941584])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    healpix_ids = _normalize_ids(healpix_ids)

    healpix_features = []
    for healpix_id in healpix_ids:
        try:
            order, _ipix, _nside = _decode_uniq(healpix_id)
            cell_polygon = _healpix_uniq_to_polygon(healpix_id)
            cell_polygon = _apply_antimeridian(cell_polygon, fix_antimeridian)
            healpix_feature = geodesic_dggs_to_feature(
                "healpix", healpix_id, order, cell_polygon, num_edges=4
            )
            healpix_features.append(healpix_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": healpix_features}


def healpix2geojson_cli():
    """
    Command-line interface for healpix2geojson supporting multiple HEALPix cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert HEALPix UNIQ cell ID(s) to GeoJSON"
    )
    parser.add_argument(
        "healpix",
        nargs="+",
        help="Input HEALPix UNIQ cell ID(s), e.g., healpix2geojson 9941583",
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
    args = parser.parse_args()
    geojson_data = json.dumps(
        healpix2geojson(args.healpix, fix_antimeridian=args.fix_antimeridian)
    )
    print(geojson_data)


if __name__ == "__main__":
    healpix2geojson_cli()

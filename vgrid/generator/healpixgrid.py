"""
HEALPix Grid Generator Module

Generates HEALPix UNIQ DGGS grids for specified resolutions with automatic
cell generation and validation.

Key Functions:
- healpix_grid(): Main grid generation function for whole world
- healpix_grid_within_bbox(): Grid generation within bounding box
- healpixgrid(): User-facing function with multiple output formats
- healpixgrid_cli(): Command-line interface for grid generation
"""

import argparse

from shapely.geometry import box
from tqdm import tqdm
import geopandas as gpd

from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    validate_bbox,
    validate_healpix_resolution,
    convert_to_output_format,
    add_verbose_argument,
)
from vgrid.conversion.dggs2geo.healpix2geo import healpix2geo
from vgrid.conversion.dggscompact.healpixcompact import healpix_compact
from vgrid.dggs.healpix import (
    nside2npix,
    order2nside,
    orderpix2uniq,
    queryBoxInclusiveNest,
    uniq2orderpix,
)


def _empty_healpix_gdf():
    return gpd.GeoDataFrame(
        columns=[
            "healpix",
            "resolution",
            "center_lat",
            "center_lon",
            "avg_edge_len",
            "cell_area",
            "cell_perimeter",
            "geometry",
        ],
        geometry="geometry",
        crs="EPSG:4326",
    )


def _healpix_row_from_uniq(uniq_id, fix_antimeridian=None):
    cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
    order, _ipix = uniq2orderpix(int(uniq_id))
    return geodesic_dggs_to_geoseries(
        "healpix", uniq_id, order, cell_polygon, num_edges=4
    )


def healpix_grid(resolution, fix_antimeridian=None, compact=False, verbose=True):
    """Generate HEALPix UNIQ cells for the whole world at a given order."""
    resolution = validate_healpix_resolution(resolution)
    nside = order2nside(resolution)
    npix = nside2npix(nside)
    if npix > MAX_CELLS:
        raise ValueError(
            f"Resolution {resolution} will generate {npix} cells which exceeds "
            f"the limit of {MAX_CELLS}"
        )

    uniq_ids = [orderpix2uniq(resolution, ipix) for ipix in range(npix)]
    if compact:
        uniq_ids = healpix_compact(uniq_ids)

    rows = []
    for uniq_id in tqdm(uniq_ids, desc="Generating HEALPix DGGS", unit=" cells", disable=not verbose):
        rows.append(_healpix_row_from_uniq(uniq_id, fix_antimeridian))
    return gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")


def healpix_grid_within_bbox(resolution, bbox, fix_antimeridian=None, compact=False, verbose=True):
    """
    Generate HEALPix UNIQ cells intersecting a bounding box.

    Uses ``queryBoxInclusiveNest`` for candidate discovery, then keeps cells that
    intersect the bbox polygon.
    """
    resolution = validate_healpix_resolution(resolution)
    bbox = validate_bbox(bbox)
    bbox_polygon = box(*bbox)

    nside = order2nside(resolution)
    nested_ids = queryBoxInclusiveNest(
        nside, (float(bbox[0]), float(bbox[1]), float(bbox[2]), float(bbox[3]))
    )
    if not nested_ids:
        return _empty_healpix_gdf()

    uniq_ids = [orderpix2uniq(resolution, int(ipix)) for ipix in nested_ids]
    if len(uniq_ids) > MAX_CELLS:
        raise ValueError(
            f"Resolution {resolution} within bounding box {bbox} will generate "
            f"{len(uniq_ids)} cells which exceeds the limit of {MAX_CELLS}"
        )

    if compact:
        uniq_ids = healpix_compact(uniq_ids)

    rows = []
    for uniq_id in tqdm(uniq_ids, desc="Generating HEALPix DGGS", unit=" cells", disable=not verbose):
        cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
        if cell_polygon is None or cell_polygon.is_empty:
            continue
        if not cell_polygon.intersects(bbox_polygon):
            continue
        order, _ipix = uniq2orderpix(int(uniq_id))
        rows.append(
            geodesic_dggs_to_geoseries(
                "healpix", uniq_id, order, cell_polygon, num_edges=4
            )
        )

    if not rows:
        return _empty_healpix_gdf()
    return gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")


def healpix_grid_ids(resolution, compact=False):
    """Return a list of HEALPix UNIQ IDs for the whole world at a given resolution."""
    resolution = validate_healpix_resolution(resolution)
    nside = order2nside(resolution)
    npix = nside2npix(nside)
    uniq_ids = [orderpix2uniq(resolution, ipix) for ipix in range(npix)]
    if compact:
        uniq_ids = healpix_compact(uniq_ids)
    return uniq_ids


def healpix_grid_within_bbox_ids(resolution, bbox, compact=False, verbose=True):
    """Return HEALPix UNIQ IDs intersecting the given bounding box."""
    gdf = healpix_grid_within_bbox(resolution, bbox, compact=compact, verbose=verbose)
    if gdf.empty:
        return []
    return gdf["healpix"].tolist()


def healpixgrid(
    resolution,
    bbox=None,
    output_format="gpd",
    fix_antimeridian=None,
    compact=False,
    verbose=True,
):
    """
    Generate HEALPix grid for pure Python usage.

    Args:
        resolution (int): HEALPix order [0..29]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat].
            Defaults to None (whole world).
        output_format (str, optional): Output format.
        fix_antimeridian (str, optional): Antimeridian fixing method.
        compact (bool, optional): Enable HEALPix compact mode to reduce cell count.

    Returns:
        Output in the requested format (GeoDataFrame, GeoJSON, file path, etc.)
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        nside = order2nside(validate_healpix_resolution(resolution))
        num_cells = nside2npix(nside)
        if num_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {num_cells} cells which "
                f"exceeds the limit of {MAX_CELLS}"
            )
        gdf = healpix_grid(
            resolution, fix_antimeridian=fix_antimeridian, compact=compact, verbose=verbose
        )
    else:
        gdf = healpix_grid_within_bbox(
            resolution,
            bbox,
            fix_antimeridian=fix_antimeridian,
            compact=compact,
            verbose=verbose,
        )
    output_name = f"healpix_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def healpixgrid_cli():
    """CLI interface for generating HEALPix grid."""
    parser = argparse.ArgumentParser(description="Generate HEALPix DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..29]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box: min_lon min_lat max_lon max_lat (default: whole world)",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable HEALPix compact mode to reduce cell count",
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

    add_verbose_argument(parser)
    args = parser.parse_args()
    try:
        result = healpixgrid(
            args.resolution,
            args.bbox,
            args.output_format,
            fix_antimeridian=args.fix_antimeridian,
            compact=args.compact,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    healpixgrid_cli()

"""
Geohash Aggregate Module

Roll Geohash cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    geohash_agg: Group Geohash IDs by parent cell at a resolution
    geohashagg: Read cells into a GeoDataFrame and aggregate by parent
    geohashagg_cli: Command-line interface
"""

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_geohash_resolution
from vgrid.conversion.dggs2geo.geohash2geo import geohash2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)


def geohash_parent_at(geohash_id, resolution):
    """Parent of ``geohash_id`` at ``resolution``: its first ``resolution`` characters."""
    geohash_id = str(geohash_id)
    check_not_coarser(geohash_id, len(geohash_id), resolution, "Geohash")
    return geohash_id[:resolution]


def geohash_agg(geohash_ids, resolution, bags=None, verbose=True):
    """Group Geohash IDs by their parent at ``resolution``.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_geohash_resolution(resolution)
    return group_by_parent(
        geohash_ids,
        lambda cid: geohash_parent_at(cid, resolution),
        "Geohash",
        bags=bags,
        verbose=verbose,
    )


def geohashagg(
    input_data,
    resolution,
    geohash_id="geohash",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate Geohash cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` on, parent geometry and
    graticule metrics are added; otherwise the result is a DataFrame of the
    Geohash id and the aggregated column.
    """
    if not geohash_id:
        geohash_id = "geohash"
    resolution = validate_geohash_resolution(resolution)

    def row_fn(parent_id):
        return dggs_cell_row(
            "geohash", parent_id, resolution, geohash2geo(parent_id), cell_metrics=True
        )

    return run_agg(
        input_data,
        geohash_id,
        "Geohash",
        "geohash",
        lambda cid: geohash_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def geohashagg_cli():
    """Command-line interface for geohashagg."""
    parser = agg_arg_parser("Geohash", "geohashcompact")
    args = parser.parse_args()
    result = geohashagg(
        args.input,
        resolution=args.resolution,
        geohash_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

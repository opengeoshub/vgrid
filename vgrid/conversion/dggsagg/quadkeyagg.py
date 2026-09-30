"""
Quadkey Aggregate Module

Roll Quadkey cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    quadkey_agg: Group Quadkeys by parent cell at a resolution
    quadkeyagg: Read cells into a GeoDataFrame and aggregate by parent
    quadkeyagg_cli: Command-line interface
"""

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_quadkey_resolution
from vgrid.conversion.dggs2geo.quadkey2geo import quadkey2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)


def quadkey_parent_at(quadkey_id, resolution):
    """Parent of ``quadkey_id`` at ``resolution``: its first ``resolution`` digits."""
    quadkey_id = str(quadkey_id)
    check_not_coarser(quadkey_id, len(quadkey_id), resolution, "Quadkey")
    return quadkey_id[:resolution]


def quadkey_agg(quadkey_ids, resolution, bags=None, verbose=True):
    """Group Quadkeys by their parent at ``resolution``.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_quadkey_resolution(resolution)
    return group_by_parent(
        quadkey_ids,
        lambda cid: quadkey_parent_at(cid, resolution),
        "Quadkey",
        bags=bags,
        verbose=verbose,
    )


def quadkeyagg(
    input_data,
    resolution,
    quadkey_id="quadkey",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate Quadkey cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` on, parent geometry and
    graticule metrics are added; otherwise the result is a DataFrame of the
    Quadkey and the aggregated column.
    """
    if not quadkey_id:
        quadkey_id = "quadkey"
    resolution = validate_quadkey_resolution(resolution)

    def row_fn(parent_id):
        return dggs_cell_row(
            "quadkey", parent_id, resolution, quadkey2geo(parent_id), cell_metrics=True
        )

    return run_agg(
        input_data,
        quadkey_id,
        "Quadkey",
        "quadkey",
        lambda cid: quadkey_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def quadkeyagg_cli():
    """Command-line interface for quadkeyagg."""
    parser = agg_arg_parser("Quadkey", "quadkeycompact")
    args = parser.parse_args()
    result = quadkeyagg(
        args.input,
        resolution=args.resolution,
        quadkey_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

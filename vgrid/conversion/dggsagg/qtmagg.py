"""
QTM Aggregate Module

Roll QTM cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    qtm_agg: Group QTM IDs by parent cell at a resolution
    qtmagg: Read cells into a GeoDataFrame and aggregate by parent
    qtmagg_cli: Command-line interface
"""

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_qtm_resolution
from vgrid.conversion.dggs2geo.qtm2geo import qtm2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)


def qtm_parent_at(qtm_id, resolution):
    """Parent of ``qtm_id`` at ``resolution``: its first ``resolution`` digits."""
    qtm_id = str(qtm_id)
    check_not_coarser(qtm_id, len(qtm_id), resolution, "QTM")
    return qtm_id[:resolution]


def qtm_agg(qtm_ids, resolution, bags=None, verbose=True):
    """Group QTM cell IDs by their parent at ``resolution``.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_qtm_resolution(resolution)
    return group_by_parent(
        qtm_ids,
        lambda cid: qtm_parent_at(cid, resolution),
        "QTM",
        bags=bags,
        verbose=verbose,
    )


def qtmagg(
    input_data,
    resolution,
    qtm_id="qtm",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate QTM cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` off the result is a
    DataFrame of the QTM id and the aggregated column.
    """
    if not qtm_id:
        qtm_id = "qtm"
    resolution = validate_qtm_resolution(resolution)

    def row_fn(parent_id):
        cell_polygon = qtm2geo(parent_id)
        return dggs_cell_row("qtm", parent_id, resolution, cell_polygon, 3, True)

    return run_agg(
        input_data,
        qtm_id,
        "QTM",
        "qtm",
        lambda cid: qtm_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def qtmagg_cli():
    """Command-line interface for qtmagg."""
    parser = agg_arg_parser("QTM", "qtmcompact")
    args = parser.parse_args()
    result = qtmagg(
        args.input,
        resolution=args.resolution,
        qtm_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

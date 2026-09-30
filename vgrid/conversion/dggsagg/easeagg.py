"""
EASE Aggregate Module

Roll EASE-DGGS cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    ease_agg: Group EASE IDs by parent cell at a resolution
    easeagg: Read cells into a GeoDataFrame and aggregate by parent
    easeagg_cli: Command-line interface
"""

import re

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_ease_resolution
from vgrid.conversion.dggs2geo.ease2geo import ease2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)

_EASE_ID = re.compile(r"L(\d+)\.(.+)")


def ease_parent_at(ease_id, resolution):
    """Parent of ``ease_id`` at ``resolution``, or the cell itself when already there.

    An EASE ID is ``L{res}.{level-0 row/col}.{one block per finer level}``, so the
    parent keeps the first ``resolution + 1`` blocks.
    """
    ease_id = str(ease_id)
    match = _EASE_ID.match(ease_id)
    if not match:
        raise ValueError(f"Invalid EASE ID <{ease_id}>.")
    check_not_coarser(ease_id, int(match.group(1)), resolution, "EASE")
    blocks = match.group(2).split(".")
    return f"L{resolution}." + ".".join(blocks[: resolution + 1])


def ease_agg(ease_ids, resolution, bags=None, verbose=True):
    """Group EASE cell IDs by their parent at ``resolution``.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_ease_resolution(resolution)
    return group_by_parent(
        ease_ids,
        lambda cid: ease_parent_at(cid, resolution),
        "EASE",
        bags=bags,
        verbose=verbose,
    )


def easeagg(
    input_data,
    resolution,
    ease_id="ease",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate EASE cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` off the result is a
    DataFrame of the EASE id and the aggregated column.
    """
    if not ease_id:
        ease_id = "ease"
    resolution = validate_ease_resolution(resolution)

    def row_fn(parent_id):
        cell_polygon = ease2geo(parent_id)
        return dggs_cell_row("ease", parent_id, resolution, cell_polygon, 4, True)

    return run_agg(
        input_data,
        ease_id,
        "EASE",
        "ease",
        lambda cid: ease_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def easeagg_cli():
    """Command-line interface for easeagg."""
    parser = agg_arg_parser("EASE", "easecompact")
    args = parser.parse_args()
    result = easeagg(
        args.input,
        resolution=args.resolution,
        ease_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

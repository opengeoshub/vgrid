"""
Tilecode Aggregate Module

Roll Tilecode cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    tilecode_agg: Group Tilecode IDs by parent cell at a resolution
    tilecodeagg: Read cells into a GeoDataFrame and aggregate by parent
    tilecodeagg_cli: Command-line interface
"""

import re

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_tilecode_resolution
from vgrid.conversion.dggs2geo.tilecode2geo import tilecode2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)

_TILECODE_ID = re.compile(r"z(\d+)x(\d+)y(\d+)")


def tilecode_parent_at(tilecode_id, resolution):
    """Parent of ``tilecode_id`` (``z{z}x{x}y{y}``) at zoom ``resolution``."""
    tilecode_id = str(tilecode_id)
    match = _TILECODE_ID.fullmatch(tilecode_id)
    if not match:
        raise ValueError(f"Invalid Tilecode ID <{tilecode_id}>.")
    z, x, y = (int(g) for g in match.groups())
    check_not_coarser(tilecode_id, z, resolution, "Tilecode")
    shift = z - resolution
    return f"z{resolution}x{x >> shift}y{y >> shift}"


def tilecode_agg(tilecode_ids, resolution, bags=None, verbose=True):
    """Group Tilecode IDs by their parent at ``resolution``.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_tilecode_resolution(resolution)
    return group_by_parent(
        tilecode_ids,
        lambda cid: tilecode_parent_at(cid, resolution),
        "Tilecode",
        bags=bags,
        verbose=verbose,
    )


def tilecodeagg(
    input_data,
    resolution,
    tilecode_id="tilecode",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate Tilecode cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` on, parent geometry and
    graticule metrics are added; otherwise the result is a DataFrame of the
    Tilecode id and the aggregated column.
    """
    if not tilecode_id:
        tilecode_id = "tilecode"
    resolution = validate_tilecode_resolution(resolution)

    def row_fn(parent_id):
        return dggs_cell_row(
            "tilecode",
            parent_id,
            resolution,
            tilecode2geo(parent_id),
            cell_metrics=True,
        )

    return run_agg(
        input_data,
        tilecode_id,
        "Tilecode",
        "tilecode",
        lambda cid: tilecode_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def tilecodeagg_cli():
    """Command-line interface for tilecodeagg."""
    parser = agg_arg_parser("Tilecode", "tilecodecompact")
    args = parser.parse_args()
    result = tilecodeagg(
        args.input,
        resolution=args.resolution,
        tilecode_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

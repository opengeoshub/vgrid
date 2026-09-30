"""
ISEA4T Aggregate Module

Roll ISEA4T cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    isea4t_agg: Group ISEA4T IDs by parent cell at a resolution
    isea4tagg: Read cells into a GeoDataFrame and aggregate by parent
    isea4tagg_cli: Command-line interface

Note: Geometry (``cell_metrics=True``) needs OpenEaggr, which is Windows-only.
"""

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_isea4t_resolution
from vgrid.utils.constants import FIX_ANTIMERIDIAN_CHOICES
from vgrid.conversion.dggs2geo.isea4t2geo import isea4t2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)


def isea4t_parent_at(isea4t_id, resolution):
    """Parent of ``isea4t_id`` at ``resolution``, or the cell itself when already there."""
    isea4t_id = str(isea4t_id)
    check_not_coarser(isea4t_id, len(isea4t_id) - 2, resolution, "ISEA4T")
    return isea4t_id[: resolution + 2]


def isea4t_agg(isea4t_ids, resolution, bags=None, verbose=True):
    """Group ISEA4T cell IDs by their parent at ``resolution``.

    Every cell is assigned to its ancestor at ``resolution``, whether or not its
    siblings are present. Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_isea4t_resolution(resolution)
    return group_by_parent(
        isea4t_ids,
        lambda cid: isea4t_parent_at(cid, resolution),
        "ISEA4T",
        bags=bags,
        verbose=verbose,
    )


def isea4tagg(
    input_data,
    resolution,
    isea4t_id="isea4t",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    cell_metrics=False,
    verbose=True,
):
    """Aggregate ISEA4T cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` off the result is a
    DataFrame of the ISEA4T id and the aggregated column.
    """
    if not isea4t_id:
        isea4t_id = "isea4t"
    resolution = validate_isea4t_resolution(resolution)

    def row_fn(parent_id):
        cell_polygon = isea4t2geo(parent_id, fix_antimeridian=fix_antimeridian)
        return dggs_cell_row("isea4t", parent_id, resolution, cell_polygon, 3, True)

    return run_agg(
        input_data,
        isea4t_id,
        "ISEA4T",
        "isea4t",
        lambda cid: isea4t_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def isea4tagg_cli():
    """Command-line interface for isea4tagg."""
    parser = agg_arg_parser("ISEA4T", "isea4tcompact")
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none.",
    )
    args = parser.parse_args()
    result = isea4tagg(
        args.input,
        resolution=args.resolution,
        isea4t_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

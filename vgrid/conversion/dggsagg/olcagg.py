"""
OLC Aggregate Module

Roll OLC (Open Location Code) cells up to a parent resolution and aggregate
values there. Incomplete child sets are included. This is a group-by parent,
not compaction.

Key Functions:
    olc_agg: Group OLC codes by parent cell at a resolution
    olcagg: Read cells into a GeoDataFrame and aggregate by parent
    olcagg_cli: Command-line interface
"""

from vgrid.dggs import olc
from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_olc_resolution
from vgrid.conversion.dggs2geo.olc2geo import olc2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    climb_to_parent,
    group_by_parent,
    print_structured,
    run_agg,
)


def get_olc_resolution(olc_id):
    return olc.decode(olc_id).codeLength


def olc_parent_at(olc_id, resolution):
    """Parent of ``olc_id`` at ``resolution`` (a valid OLC code length)."""
    return climb_to_parent(
        str(olc_id), resolution, get_olc_resolution, olc.olc_parent, "OLC"
    )


def olc_agg(olc_ids, resolution, bags=None, verbose=True):
    """Group OLC codes by their parent at ``resolution``.

    ``resolution`` must be a valid code length (2, 4, 6, 8, 10, 11-15). Cells
    coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_olc_resolution(resolution)
    return group_by_parent(
        olc_ids,
        lambda cid: olc_parent_at(cid, resolution),
        "OLC",
        bags=bags,
        verbose=verbose,
    )


def olcagg(
    input_data,
    resolution,
    olc_id="olc",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate OLC cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` on, parent geometry and
    graticule metrics are added; otherwise the result is a DataFrame of the
    OLC code and the aggregated column.
    """
    if not olc_id:
        olc_id = "olc"
    resolution = validate_olc_resolution(resolution)

    def row_fn(parent_id):
        return dggs_cell_row(
            "olc", parent_id, resolution, olc2geo(parent_id), cell_metrics=True
        )

    return run_agg(
        input_data,
        olc_id,
        "OLC",
        "olc",
        lambda cid: olc_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def olcagg_cli():
    """Command-line interface for olcagg."""
    parser = agg_arg_parser("OLC", "olccompact")
    args = parser.parse_args()
    result = olcagg(
        args.input,
        resolution=args.resolution,
        olc_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

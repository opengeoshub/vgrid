"""
DIGIPIN Aggregate Module

Roll DIGIPIN cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    digipin_agg: Group DIGIPIN codes by parent cell at a resolution
    digipinagg: Read cells into a GeoDataFrame and aggregate by parent
    digipinagg_cli: Command-line interface
"""

from vgrid.dggs.digipin import digipin_parent, digipin_resolution
from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_digipin_resolution
from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    climb_to_parent,
    group_by_parent,
    print_structured,
    run_agg,
)


def _digipin_resolution(digipin_id):
    res = digipin_resolution(digipin_id)
    if isinstance(res, str):
        raise ValueError(f"Invalid DIGIPIN <{digipin_id}>.")
    return res


def _digipin_parent(digipin_id):
    parent = digipin_parent(digipin_id)
    return None if parent == "Invalid DIGIPIN" else parent


def digipin_parent_at(digipin_id, resolution):
    """Parent of ``digipin_id`` at ``resolution``."""
    return climb_to_parent(
        str(digipin_id), resolution, _digipin_resolution, _digipin_parent, "DIGIPIN"
    )


def digipin_agg(digipin_ids, resolution, bags=None, verbose=True):
    """Group DIGIPIN codes by their parent at ``resolution``.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_digipin_resolution(resolution)
    return group_by_parent(
        digipin_ids,
        lambda cid: digipin_parent_at(cid, resolution),
        "DIGIPIN",
        bags=bags,
        verbose=verbose,
    )


def digipinagg(
    input_data,
    resolution,
    digipin_id="digipin",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate DIGIPIN cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` on, parent geometry and
    graticule metrics are added; otherwise the result is a DataFrame of the
    DIGIPIN code and the aggregated column.
    """
    if not digipin_id:
        digipin_id = "digipin"
    resolution = validate_digipin_resolution(resolution)

    def row_fn(parent_id):
        return dggs_cell_row(
            "digipin", parent_id, resolution, digipin2geo(parent_id), cell_metrics=True
        )

    return run_agg(
        input_data,
        digipin_id,
        "DIGIPIN",
        "digipin",
        lambda cid: digipin_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def digipinagg_cli():
    """Command-line interface for digipinagg."""
    parser = agg_arg_parser("DIGIPIN", "digipincompact")
    args = parser.parse_args()
    result = digipinagg(
        args.input,
        resolution=args.resolution,
        digipin_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

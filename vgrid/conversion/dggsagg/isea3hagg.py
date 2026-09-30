"""
ISEA3H Aggregate Module

Roll ISEA3H cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

ISEA3H is aperture 3, so a cell can overlap several coarser cells. Each cell is
assigned to the single parent that contains its centroid.

Key Functions:
    isea3h_agg: Group ISEA3H IDs by parent cell at a resolution
    isea3hagg: Read cells into a GeoDataFrame and aggregate by parent
    isea3hagg_cli: Command-line interface

Note: This module is only supported on Windows systems due to OpenEaggr dependency.
"""

import platform

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.dggs.eaggr.enums.model import Model

    isea3h_dggs = Eaggr(Model.ISEA3H)

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_isea3h_resolution
from vgrid.utils.constants import (
    FIX_ANTIMERIDIAN_CHOICES,
    ISEA3H_ACCURACY_RES_DICT,
    ISEA3H_RES_ACCURACY_DICT,
)
from vgrid.conversion.dggs2geo.isea3h2geo import isea3h2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)


def isea3h_parent_at(isea3h_id, resolution):
    """Cell at ``resolution`` containing the centroid of ``isea3h_id``."""
    isea3h_id = str(isea3h_id)
    centroid = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(isea3h_id))
    cell_res = ISEA3H_ACCURACY_RES_DICT.get(centroid.get_accuracy())
    if cell_res is None:
        raise ValueError(f"Invalid ISEA3H ID <{isea3h_id}>.")
    check_not_coarser(isea3h_id, cell_res, resolution, "ISEA3H")
    if cell_res == resolution:
        return isea3h_id
    point = LatLongPoint(
        centroid.get_latitude(),
        centroid.get_longitude(),
        ISEA3H_RES_ACCURACY_DICT[resolution],
    )
    return str(isea3h_dggs.convert_point_to_dggs_cell(point).get_cell_id())


def isea3h_agg(isea3h_ids, resolution, bags=None, verbose=True):
    """Group ISEA3H cell IDs by the parent at ``resolution`` containing each centroid.

    Cells coarser than ``resolution`` raise ``ValueError``.
    """
    resolution = validate_isea3h_resolution(resolution)
    return group_by_parent(
        isea3h_ids,
        lambda cid: isea3h_parent_at(cid, resolution),
        "ISEA3H",
        bags=bags,
        verbose=verbose,
    )


def isea3hagg(
    input_data,
    resolution,
    isea3h_id="isea3h",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    cell_metrics=False,
    verbose=True,
):
    """Aggregate ISEA3H cells by parent at ``resolution``.

    Parameters match ``h3agg``. With ``cell_metrics`` off the result is a
    DataFrame of the ISEA3H id and the aggregated column.
    """
    if not isea3h_id:
        isea3h_id = "isea3h"
    resolution = validate_isea3h_resolution(resolution)
    num_edges = 3 if resolution == 0 else 6

    def row_fn(parent_id):
        cell_polygon = isea3h2geo(parent_id, fix_antimeridian=fix_antimeridian)
        return dggs_cell_row(
            "isea3h", parent_id, resolution, cell_polygon, num_edges, True
        )

    return run_agg(
        input_data,
        isea3h_id,
        "ISEA3H",
        "isea3h",
        lambda cid: isea3h_parent_at(cid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def isea3hagg_cli():
    """Command-line interface for isea3hagg."""
    parser = agg_arg_parser("ISEA3H", "isea3hcompact")
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none.",
    )
    args = parser.parse_args()
    result = isea3hagg(
        args.input,
        resolution=args.resolution,
        isea3h_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

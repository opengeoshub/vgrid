"""
DGGAL Aggregate Module

Roll DGGAL zones up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Each zone is assigned to the zone at the target level that contains its
centroid. For nested DGGRS (e.g. ISEA4R, ISEA9R) that is its ancestor; for
non-nested ones (e.g. ISEA3H, ISEA7H) it picks a single parent.

Key Functions:
    dggal_agg: Group DGGAL zone IDs by parent zone at a resolution
    dggalagg: Read zones into a GeoDataFrame and aggregate by parent
    dggalagg_cli: Command-line interface
"""

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import validate_dggal_resolution, validate_dggal_type
from vgrid.utils.constants import DGGAL_TYPES
from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)

from dggal import *

app = Application(appGlobals=globals())
pydggal_setup(app)


def _dggal_dggrs(dggs_type):
    return globals()[DGGAL_TYPES[dggs_type]["class_name"]]()


def dggal_parent_at(dggrs, zone_id, resolution):
    """Zone at ``resolution`` containing the centroid of ``zone_id``."""
    zone_id = str(zone_id)
    zone = dggrs.getZoneFromTextID(zone_id)
    check_not_coarser(zone_id, dggrs.getZoneLevel(zone), resolution, "DGGAL")
    if dggrs.getZoneLevel(zone) == resolution:
        return zone_id
    parent = dggrs.getZoneFromWGS84Centroid(
        resolution, dggrs.getZoneWGS84Centroid(zone)
    )
    return dggrs.getZoneTextID(parent)


def dggal_agg(dggs_type, zone_ids, resolution, bags=None, verbose=True):
    """Group DGGAL zone IDs by the zone at ``resolution`` containing each centroid.

    Zones coarser than ``resolution`` raise ``ValueError``.
    """
    dggs_type = validate_dggal_type(dggs_type)
    resolution = validate_dggal_resolution(dggs_type, resolution)
    dggrs = _dggal_dggrs(dggs_type)
    return group_by_parent(
        zone_ids,
        lambda zid: dggal_parent_at(dggrs, zid, resolution),
        "DGGAL",
        bags=bags,
        verbose=verbose,
    )


def dggalagg(
    dggs_type,
    input_data,
    resolution,
    zone_id=None,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    split_antimeridian=False,
    cell_metrics=False,
    verbose=True,
):
    """Aggregate DGGAL zones by parent at ``resolution``.

    ``dggs_type`` is a DGGAL type such as ``isea4r`` or ``isea7h``. The zone ID
    column defaults to ``dggal_<dggs_type>``. Other parameters match
    ``h3agg``; with ``cell_metrics`` off the result is a DataFrame of the zone
    id and the aggregated column.
    """
    dggs_type = validate_dggal_type(dggs_type)
    resolution = validate_dggal_resolution(dggs_type, resolution)
    id_col = f"dggal_{dggs_type}"
    if not zone_id:
        zone_id = id_col
    dggrs = _dggal_dggrs(dggs_type)

    def row_fn(parent_id):
        zone = dggrs.getZoneFromTextID(parent_id)
        cell_polygon = dggal2geo(
            dggs_type, parent_id, split_antimeridian=split_antimeridian
        )
        return dggs_cell_row(
            id_col,
            parent_id,
            resolution,
            cell_polygon,
            dggrs.countZoneEdges(zone),
            True,
        )

    return run_agg(
        input_data,
        zone_id,
        "DGGAL",
        f"dggal_{dggs_type}",
        lambda zid: dggal_parent_at(dggrs, zid, resolution),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def dggalagg_cli():
    """Command-line interface for dggalagg."""
    parser = agg_arg_parser("DGGAL", "dggalcompact")
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL type",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable Antimeridian splitting",
    )
    args = parser.parse_args()
    result = dggalagg(
        args.dggs_type,
        args.input,
        resolution=args.resolution,
        zone_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        split_antimeridian=args.split_antimeridian,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

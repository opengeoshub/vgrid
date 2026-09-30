"""
rHEALPix Aggregate Module

Roll rHEALPix cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    rhealpix_agg: Group rHEALPix IDs by parent cell at a resolution
    rhealpixagg: Read cells into a GeoDataFrame and aggregate by parent
    rhealpixagg_cli: Command-line interface
"""

from vgrid.utils.geometry import dggs_cell_row
from vgrid.utils.io import (
    add_rhealpix_n_side_argument,
    get_rhealpix_dggs,
    rhealpix_cell_from_id,
    validate_rhealpix_resolution,
)
from vgrid.utils.constants import FIX_ANTIMERIDIAN_CHOICES
from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    check_not_coarser,
    group_by_parent,
    print_structured,
    run_agg,
)


def rhealpix_parent_at(rhealpix_id, resolution, dggs):
    """Parent of ``rhealpix_id`` at ``resolution``, or the cell itself when already there."""
    suid = dggs.parse_index(str(rhealpix_id))
    if suid is None:
        raise ValueError(f"Invalid rHEALPix ID <{rhealpix_id}>.")
    check_not_coarser(rhealpix_id, len(suid) - 1, resolution, "rHEALPix")
    return dggs.format_index(suid[: resolution + 1])


def rhealpix_agg(rhealpix_ids, resolution, bags=None, verbose=True, N_side=3):
    """Group rHEALPix cell IDs by their parent at ``resolution``.

    Every cell is assigned to its ancestor at ``resolution``, whether or not its
    siblings are present. Cells already at ``resolution`` stay put. Cells
    coarser than ``resolution`` raise ``ValueError``. ``bags`` works as in
    ``h3_agg``.
    """
    resolution = validate_rhealpix_resolution(resolution)
    dggs = get_rhealpix_dggs(N_side=N_side)
    return group_by_parent(
        rhealpix_ids,
        lambda rid: rhealpix_parent_at(rid, resolution, dggs),
        "rHEALPix",
        bags=bags,
        verbose=verbose,
    )


def rhealpixagg(
    input_data,
    resolution,
    rhealpix_id="rhealpix",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    cell_metrics=False,
    verbose=True,
    N_side=3,
):
    """Aggregate rHEALPix cells by parent at ``resolution``.

    Parameters match ``h3agg``; ``N_side`` (2 or 3) selects the rHEALPix grid.
    With ``cell_metrics`` off the result is a DataFrame of the rHEALPix id and
    the aggregated column; with it on, parent geometry and geodesic metrics
    are added.
    """
    if not rhealpix_id:
        rhealpix_id = "rhealpix"
    resolution = validate_rhealpix_resolution(resolution)
    dggs = get_rhealpix_dggs(N_side=N_side)

    def row_fn(parent_id):
        cell_polygon = rhealpix2geo(
            parent_id, fix_antimeridian=fix_antimeridian, N_side=N_side
        )
        cell = rhealpix_cell_from_id(parent_id, dggs=dggs)
        num_edges = 3 if cell.ellipsoidal_shape == "dart" else 4
        return dggs_cell_row(
            "rhealpix", parent_id, resolution, cell_polygon, num_edges, True
        )

    return run_agg(
        input_data,
        rhealpix_id,
        "rHEALPix",
        "rhealpix",
        lambda rid: rhealpix_parent_at(rid, resolution, dggs),
        row_fn,
        agg=agg,
        numeric_col=numeric_col,
        output_format=output_format,
        cell_metrics=cell_metrics,
        verbose=verbose,
    )


def rhealpixagg_cli():
    """Command-line interface for rhealpixagg."""
    parser = agg_arg_parser("rHEALPix", "rhealpixcompact")
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none.",
    )
    add_rhealpix_n_side_argument(parser)
    args = parser.parse_args()
    result = rhealpixagg(
        args.input,
        resolution=args.resolution,
        rhealpix_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
        N_side=args.N_side,
    )
    print_structured(args.output_format, result)

"""
ISEA4T vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the ISEA4T layer in Vgrid Viz. See ``dggs2cogp.common``.
Parent polygons need OpenEaggr, which is Windows-only.

Key Functions:
    isea4t2cogp: Convert an ISEA4T vector layer into Cloud Optimized GeoParquet
    isea4t2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.isea4t2geo import isea4t2geo
from vgrid.conversion.dggsagg.isea4tagg import isea4t_parent_at
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    add_fix_antimeridian_argument,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    normalize_fix_antimeridian,
    per_cell_geometries,
    per_cell_parents,
    run_cli,
    scale_for_zoom,
)
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import get_isea4t_resolution_from_scale_denominator

DEFAULT_ISEA4T_COL = "isea4t"
DEFAULT_MIN_RES = DGGS_TYPES["isea4t"]["min_res"]
DEFAULT_FIX_ANTIMERIDIAN = "shift"


def _resolution_for_zoom(zoom):
    return get_isea4t_resolution_from_scale_denominator(
        scale_for_zoom(zoom), relative_depth=8, mm_per_pixel=0.28
    )


def isea4t_spec(fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN):
    return CogpSpec(
        label="ISEA4T",
        min_res=DGGS_TYPES["isea4t"]["min_res"],
        max_res=DGGS_TYPES["isea4t"]["max_res"],
        cell_resolution=lambda c: len(c) - 2,
        parents=per_cell_parents(isea4t_parent_at),
        geometries=per_cell_geometries(
            lambda c, _res: isea4t2geo(c, fix_antimeridian=fix_antimeridian)
        ),
        resolution_for_zoom=_resolution_for_zoom,
    )


def isea4t2cogp(
    input_path,
    output_path,
    agg_col=None,
    isea4t_col=DEFAULT_ISEA4T_COL,
    min_res=DEFAULT_MIN_RES,
    fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN,
    verbose=True,
):
    """Convert an ISEA4T vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        isea4t_spec(normalize_fix_antimeridian(fix_antimeridian)),
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=isea4t_col,
        min_res=min_res,
        verbose=verbose,
    )


def isea4t2cogp_cli():
    parser = cogp_arg_parser(
        "ISEA4T", "isea4t_col", DEFAULT_ISEA4T_COL, DEFAULT_MIN_RES, "isea4t_12"
    )
    add_fix_antimeridian_argument(parser, DEFAULT_FIX_ANTIMERIDIAN)
    args = parser.parse_args()
    run_cli(
        lambda: isea4t2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            fix_antimeridian=args.fix_antimeridian,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    isea4t2cogp_cli()

"""
rHEALPix vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the rHEALPix layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    rhealpix2cogp: Convert an rHEALPix vector layer into Cloud Optimized GeoParquet
    rhealpix2cogp_cli: Command-line interface
"""

from functools import lru_cache

from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo
from vgrid.conversion.dggsagg.rhealpixagg import rhealpix_parent_at
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    add_fix_antimeridian_argument,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    normalize_fix_antimeridian,
    per_cell_geometries,
    run_cli,
    scale_for_zoom,
)
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import get_rhealpix_resolution_from_scale_denominator
from vgrid.utils.io import add_rhealpix_n_side_argument, get_rhealpix_dggs

DEFAULT_RHEALPIX_COL = "rhealpix"
DEFAULT_MIN_RES = DGGS_TYPES["rhealpix"]["min_res"]
DEFAULT_FIX_ANTIMERIDIAN = "shift"


@lru_cache(maxsize=None)
def _resolution_for_zoom_fn(N_side):
    relative_depth = 8 if N_side == 2 else 5

    def resolution_for_zoom(zoom):
        return get_rhealpix_resolution_from_scale_denominator(
            scale_for_zoom(zoom),
            relative_depth=relative_depth,
            mm_per_pixel=0.28,
            N_side=N_side,
        )

    return resolution_for_zoom


def rhealpix_spec(N_side=3, fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN):
    dggs = get_rhealpix_dggs(N_side=N_side)
    return CogpSpec(
        label="rHEALPix",
        min_res=DGGS_TYPES["rhealpix"]["min_res"],
        max_res=DGGS_TYPES["rhealpix"]["max_res"],
        cell_resolution=lambda c: len(dggs.parse_index(c)) - 1,
        parents=lambda cells, res: [rhealpix_parent_at(c, res, dggs) for c in cells],
        geometries=per_cell_geometries(
            lambda c, _res: rhealpix2geo(
                c, fix_antimeridian=fix_antimeridian, N_side=N_side
            )
        ),
        resolution_for_zoom=_resolution_for_zoom_fn(N_side),
    )


def rhealpix2cogp(
    input_path,
    output_path,
    agg_col=None,
    rhealpix_col=DEFAULT_RHEALPIX_COL,
    min_res=DEFAULT_MIN_RES,
    fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN,
    N_side=3,
    verbose=True,
):
    """Convert an rHEALPix vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``; ``N_side`` (2 or 3) selects the rHEALPix grid.
    """
    dggs2cogp(
        rhealpix_spec(N_side, normalize_fix_antimeridian(fix_antimeridian)),
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=rhealpix_col,
        min_res=min_res,
        verbose=verbose,
    )


def rhealpix2cogp_cli():
    parser = cogp_arg_parser(
        "rHEALPix", "rhealpix_col", DEFAULT_RHEALPIX_COL, DEFAULT_MIN_RES, "rhealpix_8"
    )
    add_fix_antimeridian_argument(parser, DEFAULT_FIX_ANTIMERIDIAN)
    add_rhealpix_n_side_argument(parser)
    args = parser.parse_args()
    run_cli(
        lambda: rhealpix2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            fix_antimeridian=args.fix_antimeridian,
            N_side=args.N_side,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    rhealpix2cogp_cli()

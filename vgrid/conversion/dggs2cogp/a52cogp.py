"""
A5 vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the A5 layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    a52cogp: Convert an A5 vector layer into Cloud Optimized GeoParquet
    a52cogp_cli: Command-line interface
"""

import a5

from vgrid.conversion.dggs2geo.a52geo import a52geo
from vgrid.conversion.dggsagg.a5agg import a5_parent_at
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    add_split_antimeridian_argument,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    per_cell_geometries,
    per_cell_parents,
    run_cli,
    scale_for_zoom,
)
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import get_a5_resolution_from_scale_denominator

DEFAULT_A5_COL = "a5"
DEFAULT_MIN_RES = DGGS_TYPES["a5"]["min_res"]


def _resolution_for_zoom(zoom):
    return get_a5_resolution_from_scale_denominator(
        scale_for_zoom(zoom), relative_depth=8, mm_per_pixel=0.28
    )


def a5_spec(split_antimeridian=False):
    return CogpSpec(
        label="A5",
        min_res=DGGS_TYPES["a5"]["min_res"],
        max_res=DGGS_TYPES["a5"]["max_res"],
        cell_resolution=lambda c: a5.get_resolution(a5.hex_to_u64(c)),
        parents=per_cell_parents(a5_parent_at),
        geometries=per_cell_geometries(
            lambda c, _res: a52geo(c, split_antimeridian=split_antimeridian)
        ),
        resolution_for_zoom=_resolution_for_zoom,
    )


def a52cogp(
    input_path,
    output_path,
    agg_col=None,
    a5_col=DEFAULT_A5_COL,
    min_res=DEFAULT_MIN_RES,
    split_antimeridian=False,
    verbose=True,
):
    """Convert an A5 vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``; A5 parents use ``split_antimeridian``
    instead of ``fix_antimeridian``.
    """
    dggs2cogp(
        a5_spec(split_antimeridian),
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=a5_col,
        min_res=min_res,
        verbose=verbose,
    )


def a52cogp_cli():
    parser = cogp_arg_parser("A5", "a5_col", DEFAULT_A5_COL, DEFAULT_MIN_RES, "a5_10")
    add_split_antimeridian_argument(parser)
    args = parser.parse_args()
    run_cli(
        lambda: a52cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            split_antimeridian=args.split_antimeridian,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    a52cogp_cli()

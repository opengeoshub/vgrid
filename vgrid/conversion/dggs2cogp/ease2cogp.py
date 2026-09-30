"""
EASE-DGGS vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the EASE layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    ease2cogp: Convert an EASE vector layer into Cloud Optimized GeoParquet
    ease2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.ease2geo import ease2geo
from vgrid.conversion.dggsagg.easeagg import ease_parent_at
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    per_cell_geometries,
    per_cell_parents,
    run_cli,
)
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import get_ease_resolution

DEFAULT_EASE_COL = "ease"
DEFAULT_MIN_RES = DGGS_TYPES["ease"]["min_res"]
# (first zoom, resolution) steps from the EASE layer in Vgrid Viz.
_ZOOM_STEPS = ((23, 6), (20, 5), (16, 4), (14, 3), (12, 2), (10, 1))


def _resolution_for_zoom(zoom):
    for first_zoom, resolution in _ZOOM_STEPS:
        if zoom >= first_zoom:
            return resolution
    return 0


EASE_SPEC = CogpSpec(
    label="EASE",
    min_res=DGGS_TYPES["ease"]["min_res"],
    max_res=DGGS_TYPES["ease"]["max_res"],
    cell_resolution=get_ease_resolution,
    parents=per_cell_parents(ease_parent_at),
    geometries=per_cell_geometries(lambda c, _res: ease2geo(c)),
    resolution_for_zoom=_resolution_for_zoom,
)


def ease2cogp(
    input_path,
    output_path,
    agg_col=None,
    ease_col=DEFAULT_EASE_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert an EASE vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        EASE_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=ease_col,
        min_res=min_res,
        verbose=verbose,
    )


def ease2cogp_cli():
    parser = cogp_arg_parser("EASE", "ease_col", DEFAULT_EASE_COL, DEFAULT_MIN_RES, "ease_4")
    args = parser.parse_args()
    run_cli(
        lambda: ease2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    ease2cogp_cli()

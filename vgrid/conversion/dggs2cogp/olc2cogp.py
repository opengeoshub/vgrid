"""
OLC (Open Location Code) vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last. Only
valid code lengths (2, 4, 6, 8, 10, 11-15) become levels. Level zooms follow
the OLC layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    olc2cogp: Convert an OLC vector layer into Cloud Optimized GeoParquet
    olc2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.olc2geo import olc2geo
from vgrid.conversion.dggsagg.olcagg import get_olc_resolution, olc_parent_at
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

DEFAULT_OLC_COL = "olc"
DEFAULT_MIN_RES = DGGS_TYPES["olc"]["min_res"]
OLC_RESOLUTIONS = [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]
# (last zoom, code length) steps from the OLC layer in Vgrid Viz.
_ZOOM_STEPS = ((6, 2), (10, 4), (14, 6), (18, 8), (21, 10), (23, 11), (25, 12), (27, 13), (29, 14))


def _resolution_for_zoom(zoom):
    for last_zoom, resolution in _ZOOM_STEPS:
        if zoom <= last_zoom:
            return resolution
    return 15


OLC_SPEC = CogpSpec(
    label="OLC",
    min_res=DGGS_TYPES["olc"]["min_res"],
    max_res=DGGS_TYPES["olc"]["max_res"],
    cell_resolution=get_olc_resolution,
    parents=per_cell_parents(olc_parent_at),
    geometries=per_cell_geometries(lambda c, _res: olc2geo(c)),
    resolution_for_zoom=_resolution_for_zoom,
    valid_resolutions=OLC_RESOLUTIONS,
)


def olc2cogp(
    input_path,
    output_path,
    agg_col=None,
    olc_col=DEFAULT_OLC_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert an OLC vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        OLC_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=olc_col,
        min_res=min_res,
        verbose=verbose,
    )


def olc2cogp_cli():
    parser = cogp_arg_parser("OLC", "olc_col", DEFAULT_OLC_COL, DEFAULT_MIN_RES, "olc_10")
    args = parser.parse_args()
    run_cli(
        lambda: olc2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    olc2cogp_cli()

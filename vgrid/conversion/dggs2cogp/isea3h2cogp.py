"""
ISEA3H vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last. ISEA3H
is aperture 3, so each cell is summed into the parent containing its centroid.
Level zooms follow the ISEA3H layer in Vgrid Viz. See ``dggs2cogp.common``.
Only supported on Windows (OpenEaggr).

Key Functions:
    isea3h2cogp: Convert an ISEA3H vector layer into Cloud Optimized GeoParquet
    isea3h2cogp_cli: Command-line interface
"""

import importlib

from vgrid.conversion.dggs2geo.isea3h2geo import isea3h2geo
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
from vgrid.utils.constants import DGGS_TYPES, ISEA3H_ACCURACY_RES_DICT
from vgrid.utils.geometry import get_isea3h_resolution_from_scale_denominator

# The package re-exports the ``isea3hagg`` function under the module's name.
_isea3hagg = importlib.import_module("vgrid.conversion.dggsagg.isea3hagg")

DEFAULT_ISEA3H_COL = "isea3h"
DEFAULT_MIN_RES = DGGS_TYPES["isea3h"]["min_res"]
DEFAULT_FIX_ANTIMERIDIAN = "shift"


def _resolution_for_zoom(zoom):
    return get_isea3h_resolution_from_scale_denominator(
        scale_for_zoom(zoom), relative_depth=10, mm_per_pixel=0.28
    )


def _cell_resolution(isea3h_id):
    point = _isea3hagg.isea3h_dggs.convert_dggs_cell_to_point(
        _isea3hagg.DggsCell(isea3h_id)
    )
    return ISEA3H_ACCURACY_RES_DICT.get(point.get_accuracy())


def isea3h_spec(fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN):
    return CogpSpec(
        label="ISEA3H",
        min_res=DGGS_TYPES["isea3h"]["min_res"],
        max_res=DGGS_TYPES["isea3h"]["max_res"],
        cell_resolution=_cell_resolution,
        parents=per_cell_parents(_isea3hagg.isea3h_parent_at),
        geometries=per_cell_geometries(
            lambda c, _res: isea3h2geo(c, fix_antimeridian=fix_antimeridian)
        ),
        resolution_for_zoom=_resolution_for_zoom,
    )


def isea3h2cogp(
    input_path,
    output_path,
    agg_col=None,
    isea3h_col=DEFAULT_ISEA3H_COL,
    min_res=DEFAULT_MIN_RES,
    fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN,
    verbose=True,
):
    """Convert an ISEA3H vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        isea3h_spec(normalize_fix_antimeridian(fix_antimeridian)),
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=isea3h_col,
        min_res=min_res,
        verbose=verbose,
    )


def isea3h2cogp_cli():
    parser = cogp_arg_parser(
        "ISEA3H", "isea3h_col", DEFAULT_ISEA3H_COL, DEFAULT_MIN_RES, "isea3h_14"
    )
    add_fix_antimeridian_argument(parser, DEFAULT_FIX_ANTIMERIDIAN)
    args = parser.parse_args()
    run_cli(
        lambda: isea3h2cogp(
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
    isea3h2cogp_cli()

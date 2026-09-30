"""
DGGAL vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent zones with -agg_col summed from
their children, down to -min_res; the file resolution is written last. Each
zone is summed into the parent zone containing its centroid. Level zooms
follow the matching DGGAL layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    dggal2cogp: Convert a DGGAL vector layer into Cloud Optimized GeoParquet
    dggal2cogp_cli: Command-line interface
"""

import importlib
from functools import lru_cache

from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    add_split_antimeridian_argument,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    per_cell_geometries,
    run_cli,
    scale_for_zoom,
)
from vgrid.utils.constants import DGGAL_TYPES
from vgrid.utils.io import validate_dggal_type

# The package re-exports the ``dggalagg`` function under the module's name.
_dggalagg = importlib.import_module("vgrid.conversion.dggsagg.dggalagg")

DEFAULT_MIN_RES = 0


def _relative_depth(dggs_type):
    """Relative depth the Vgrid Viz layer passes to getLevelFromScaleDenominator."""
    if dggs_type.endswith("3h"):
        return 10
    if dggs_type.endswith("9r") or dggs_type == "rhealpix":
        return 5
    if "7h" in dggs_type:
        return 6
    return 8


@lru_cache(maxsize=None)
def _resolution_for_zoom_fn(dggs_type):
    dggrs = _dggalagg._dggal_dggrs(dggs_type)
    relative_depth = _relative_depth(dggs_type)

    def resolution_for_zoom(zoom):
        return dggrs.getLevelFromScaleDenominator(
            scale_for_zoom(zoom), relativeDepth=relative_depth, mmPerPixel=0.28
        )

    return resolution_for_zoom


def dggal_spec(dggs_type, split_antimeridian=False):
    dggs_type = validate_dggal_type(dggs_type)
    dggrs = _dggalagg._dggal_dggrs(dggs_type)
    return CogpSpec(
        label=f"DGGAL {dggs_type}",
        min_res=int(DGGAL_TYPES[dggs_type]["min_res"]),
        max_res=int(DGGAL_TYPES[dggs_type]["max_res"]),
        cell_resolution=lambda c: dggrs.getZoneLevel(dggrs.getZoneFromTextID(c)),
        parents=lambda cells, res: [
            _dggalagg.dggal_parent_at(dggrs, c, res) for c in cells
        ],
        geometries=per_cell_geometries(
            lambda c, _res: dggal2geo(
                dggs_type, c, split_antimeridian=split_antimeridian
            )
        ),
        resolution_for_zoom=_resolution_for_zoom_fn(dggs_type),
    )


def dggal2cogp(
    dggs_type,
    input_path,
    output_path,
    agg_col=None,
    dggal_col=None,
    min_res=DEFAULT_MIN_RES,
    split_antimeridian=False,
    verbose=True,
):
    """Convert a DGGAL vector layer into Cloud Optimized GeoParquet.

    ``dggs_type`` is a DGGAL type such as ``isea4r`` or ``isea7h``; the zone
    column defaults to ``dggal_<dggs_type>``. Other parameters match ``h32cogp``.
    """
    dggs_type = validate_dggal_type(dggs_type)
    dggs2cogp(
        dggal_spec(dggs_type, split_antimeridian),
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=dggal_col or f"dggal_{dggs_type}",
        min_res=min_res,
        verbose=verbose,
    )


def dggal2cogp_cli():
    parser = cogp_arg_parser(
        "DGGAL", "dggal_col", None, DEFAULT_MIN_RES, "dggal_isea4r_10"
    )
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL type",
    )
    add_split_antimeridian_argument(parser)
    args = parser.parse_args()
    run_cli(
        lambda: dggal2cogp(
            args.dggs_type,
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
    dggal2cogp_cli()

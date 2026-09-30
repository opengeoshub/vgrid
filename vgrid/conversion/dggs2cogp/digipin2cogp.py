"""
DIGIPIN vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the DIGIPIN layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    digipin2cogp: Convert a DIGIPIN vector layer into Cloud Optimized GeoParquet
    digipin2cogp_cli: Command-line interface
"""

import importlib

from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo
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

# The package re-exports the ``digipinagg`` function under the module's name.
_digipinagg = importlib.import_module("vgrid.conversion.dggsagg.digipinagg")

DEFAULT_DIGIPIN_COL = "digipin"
DEFAULT_MIN_RES = DGGS_TYPES["digipin"]["min_res"]


def _resolution_for_zoom(zoom):
    """DIGIPIN length 1 below zoom 5, then one more character every 2 zooms, up to 10."""
    if zoom < 5:
        return 1
    return min(10, 2 + int((zoom - 5) // 2))


DIGIPIN_SPEC = CogpSpec(
    label="DIGIPIN",
    min_res=DGGS_TYPES["digipin"]["min_res"],
    max_res=DGGS_TYPES["digipin"]["max_res"],
    cell_resolution=_digipinagg._digipin_resolution,
    parents=per_cell_parents(_digipinagg.digipin_parent_at),
    geometries=per_cell_geometries(lambda c, _res: digipin2geo(c)),
    resolution_for_zoom=_resolution_for_zoom,
)


def digipin2cogp(
    input_path,
    output_path,
    agg_col=None,
    digipin_col=DEFAULT_DIGIPIN_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert a DIGIPIN vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        DIGIPIN_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=digipin_col,
        min_res=min_res,
        verbose=verbose,
    )


def digipin2cogp_cli():
    parser = cogp_arg_parser(
        "DIGIPIN", "digipin_col", DEFAULT_DIGIPIN_COL, DEFAULT_MIN_RES, "digipin_8"
    )
    args = parser.parse_args()
    run_cli(
        lambda: digipin2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    digipin2cogp_cli()

"""
Quadkey vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Quadkeys are the same tiles as Tilecode, so level zooms follow the Tilecode
layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    quadkey2cogp: Convert a Quadkey vector layer into Cloud Optimized GeoParquet
    quadkey2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.quadkey2geo import quadkey2geo
from vgrid.conversion.dggsagg.quadkeyagg import quadkey_parent_at
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    per_cell_geometries,
    per_cell_parents,
    run_cli,
)
from vgrid.conversion.dggs2cogp.tilecode2cogp import tile_resolution_for_zoom
from vgrid.utils.constants import DGGS_TYPES

DEFAULT_QUADKEY_COL = "quadkey"
DEFAULT_MIN_RES = DGGS_TYPES["quadkey"]["min_res"]

QUADKEY_SPEC = CogpSpec(
    label="Quadkey",
    min_res=DGGS_TYPES["quadkey"]["min_res"],
    max_res=DGGS_TYPES["quadkey"]["max_res"],
    cell_resolution=len,
    parents=per_cell_parents(quadkey_parent_at),
    geometries=per_cell_geometries(lambda c, _res: quadkey2geo(c)),
    resolution_for_zoom=tile_resolution_for_zoom,
)


def quadkey2cogp(
    input_path,
    output_path,
    agg_col=None,
    quadkey_col=DEFAULT_QUADKEY_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert a Quadkey vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        QUADKEY_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=quadkey_col,
        min_res=min_res,
        verbose=verbose,
    )


def quadkey2cogp_cli():
    parser = cogp_arg_parser(
        "Quadkey", "quadkey_col", DEFAULT_QUADKEY_COL, DEFAULT_MIN_RES, "quadkey_14"
    )
    args = parser.parse_args()
    run_cli(
        lambda: quadkey2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    quadkey2cogp_cli()

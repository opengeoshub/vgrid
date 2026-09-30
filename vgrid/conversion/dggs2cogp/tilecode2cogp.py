"""
Tilecode vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the Tilecode layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    tilecode2cogp: Convert a Tilecode vector layer into Cloud Optimized GeoParquet
    tilecode2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.tilecode2geo import tilecode2geo
from vgrid.conversion.dggsagg.tilecodeagg import tilecode_parent_at
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    per_cell_geometries,
    per_cell_parents,
    run_cli,
    scale_for_zoom,
)
from vgrid.dggs.tilecode import tilecode_resolution
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import get_tilecode_resolution_from_scale_denominator

DEFAULT_TILECODE_COL = "tilecode"
DEFAULT_MIN_RES = DGGS_TYPES["tilecode"]["min_res"]


def tile_resolution_for_zoom(zoom):
    return get_tilecode_resolution_from_scale_denominator(
        scale_for_zoom(zoom), relative_depth=8, mm_per_pixel=0.28
    )


TILECODE_SPEC = CogpSpec(
    label="Tilecode",
    min_res=DGGS_TYPES["tilecode"]["min_res"],
    max_res=DGGS_TYPES["tilecode"]["max_res"],
    cell_resolution=tilecode_resolution,
    parents=per_cell_parents(tilecode_parent_at),
    geometries=per_cell_geometries(lambda c, _res: tilecode2geo(c)),
    resolution_for_zoom=tile_resolution_for_zoom,
)


def tilecode2cogp(
    input_path,
    output_path,
    agg_col=None,
    tilecode_col=DEFAULT_TILECODE_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert a Tilecode vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        TILECODE_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=tilecode_col,
        min_res=min_res,
        verbose=verbose,
    )


def tilecode2cogp_cli():
    parser = cogp_arg_parser(
        "Tilecode", "tilecode_col", DEFAULT_TILECODE_COL, DEFAULT_MIN_RES, "tilecode_14"
    )
    args = parser.parse_args()
    run_cli(
        lambda: tilecode2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    tilecode2cogp_cli()

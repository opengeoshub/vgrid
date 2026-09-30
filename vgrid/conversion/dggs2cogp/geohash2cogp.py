"""
Geohash vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the Geohash layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    geohash2cogp: Convert a Geohash vector layer into Cloud Optimized GeoParquet
    geohash2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.geohash2geo import geohash2geo
from vgrid.conversion.dggsagg.geohashagg import geohash_parent_at
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
from vgrid.utils.constants import DGGS_TYPES
from vgrid.utils.geometry import get_geohash_resolution_from_scale_denominator

DEFAULT_GEOHASH_COL = "geohash"
DEFAULT_MIN_RES = DGGS_TYPES["geohash"]["min_res"]


def _resolution_for_zoom(zoom):
    return get_geohash_resolution_from_scale_denominator(
        scale_for_zoom(zoom), relative_depth=3, mm_per_pixel=0.28
    )


GEOHASH_SPEC = CogpSpec(
    label="Geohash",
    min_res=DGGS_TYPES["geohash"]["min_res"],
    max_res=DGGS_TYPES["geohash"]["max_res"],
    cell_resolution=len,
    parents=per_cell_parents(geohash_parent_at),
    geometries=per_cell_geometries(lambda c, _res: geohash2geo(c)),
    resolution_for_zoom=_resolution_for_zoom,
)


def geohash2cogp(
    input_path,
    output_path,
    agg_col=None,
    geohash_col=DEFAULT_GEOHASH_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert a Geohash vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        GEOHASH_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=geohash_col,
        min_res=min_res,
        verbose=verbose,
    )


def geohash2cogp_cli():
    parser = cogp_arg_parser(
        "Geohash", "geohash_col", DEFAULT_GEOHASH_COL, DEFAULT_MIN_RES, "geohash_6"
    )
    args = parser.parse_args()
    run_cli(
        lambda: geohash2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    geohash2cogp_cli()

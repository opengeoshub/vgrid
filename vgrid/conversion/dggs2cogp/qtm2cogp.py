"""
QTM vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last.
Level zooms follow the QTM layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    qtm2cogp: Convert a QTM vector layer into Cloud Optimized GeoParquet
    qtm2cogp_cli: Command-line interface
"""

from vgrid.conversion.dggs2geo.qtm2geo import qtm2geo
from vgrid.conversion.dggsagg.qtmagg import qtm_parent_at
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
from vgrid.utils.geometry import get_qtm_resolution_from_scale_denominator

DEFAULT_QTM_COL = "qtm"
DEFAULT_MIN_RES = DGGS_TYPES["qtm"]["min_res"]


def _resolution_for_zoom(zoom):
    return get_qtm_resolution_from_scale_denominator(
        scale_for_zoom(zoom), relative_depth=8, mm_per_pixel=0.28
    )


QTM_SPEC = CogpSpec(
    label="QTM",
    min_res=DGGS_TYPES["qtm"]["min_res"],
    max_res=DGGS_TYPES["qtm"]["max_res"],
    cell_resolution=len,
    parents=per_cell_parents(qtm_parent_at),
    geometries=per_cell_geometries(lambda c, _res: qtm2geo(c)),
    resolution_for_zoom=_resolution_for_zoom,
)


def qtm2cogp(
    input_path,
    output_path,
    agg_col=None,
    qtm_col=DEFAULT_QTM_COL,
    min_res=DEFAULT_MIN_RES,
    verbose=True,
):
    """Convert a QTM vector layer into Cloud Optimized GeoParquet.

    Parameters match ``h32cogp``.
    """
    dggs2cogp(
        QTM_SPEC,
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=qtm_col,
        min_res=min_res,
        verbose=verbose,
    )


def qtm2cogp_cli():
    parser = cogp_arg_parser("QTM", "qtm_col", DEFAULT_QTM_COL, DEFAULT_MIN_RES, "qtm_12")
    args = parser.parse_args()
    run_cli(
        lambda: qtm2cogp(
            args.input,
            default_output_path(args),
            args.agg_col,
            args.id_col,
            args.min_res,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    qtm2cogp_cli()

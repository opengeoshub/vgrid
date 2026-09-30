"""
DGGRID vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells with -agg_col summed from
their children, down to -min_res; the file resolution is written last. Each
cell is summed into the parent containing its centroid. DGGRID SEQNUMs do not
carry their resolution, so the input resolution must be given. Level zooms
follow the DGGRID layer in Vgrid Viz. See ``dggs2cogp.common``.

Key Functions:
    dggrid2cogp: Convert a DGGRID vector layer into Cloud Optimized GeoParquet
    dggrid2cogp_cli: Command-line interface
"""

import importlib
import json
from functools import lru_cache
from math import floor

from vgrid.conversion.dggs2geo.dggrid2geo import dggrid2geo
from vgrid.conversion.dggs2cogp.common import (
    CogpSpec,
    add_split_antimeridian_argument,
    cogp_arg_parser,
    default_output_path,
    dggs2cogp,
    run_cli,
)
from vgrid.utils.constants import DGGRID_TYPES
from vgrid.utils.io import create_dggrid_instance, validate_dggrid_type

# The package does not export dggridagg; import the module directly.
_dggridagg = importlib.import_module("vgrid.conversion.dggsagg.dggridagg")

DEFAULT_MIN_RES = 0
# res = floor(zoom * k) per type, as in the DGGRID layer of Vgrid Viz.
_ZOOM_K = {
    "ISEA4T": 0.95,
    "FULLER4T": 0.95,
    "ISEA4D": 0.95,
    "FULLER4D": 0.95,
    "ISEA3H": 1.15,
    "FULLER3H": 1.15,
    "ISEA4H": 0.95,
    "FULLER4H": 0.95,
    "IGEO7": 0.65,
    "ISEA7H": 0.65,
    "FULLER7H": 0.65,
}


@lru_cache(maxsize=None)
def _resolution_for_zoom_fn(dggs_type):
    k = _ZOOM_K.get(dggs_type, 1.0)
    return lambda zoom: int(floor(zoom * k))


def dggrid_spec(
    dggrid_instance, dggs_type, source_resolution, split_antimeridian=False, options=None
):
    dggs_type = validate_dggrid_type(dggs_type)
    id_col = f"dggrid_{dggs_type.lower()}"

    def parents(cells, resolution):
        mapping = _dggridagg.dggrid_parent_map(
            dggrid_instance, dggs_type, cells, source_resolution, resolution
        )
        return [mapping[str(c)] for c in cells]

    def geometries(cells, resolution):
        gdf = dggrid2geo(
            dggrid_instance,
            dggs_type,
            list(cells),
            resolution,
            split_antimeridian=split_antimeridian,
            options=options,
        )
        by_id = {str(cid): geom for cid, geom in zip(gdf[id_col], gdf.geometry)}
        return [by_id.get(str(c)) for c in cells]

    return CogpSpec(
        label=f"DGGRID {dggs_type}",
        min_res=int(DGGRID_TYPES[dggs_type]["min_res"]),
        max_res=int(DGGRID_TYPES[dggs_type]["max_res"]),
        cell_resolution=None,
        parents=parents,
        geometries=geometries,
        resolution_for_zoom=_resolution_for_zoom_fn(dggs_type),
    )


def dggrid2cogp(
    dggrid_instance,
    dggs_type,
    input_path,
    output_path,
    input_resolution,
    agg_col=None,
    dggrid_col=None,
    min_res=DEFAULT_MIN_RES,
    split_antimeridian=False,
    options=None,
    verbose=True,
):
    """Convert a DGGRID vector layer at ``input_resolution`` into Cloud Optimized GeoParquet.

    The cell column defaults to ``dggrid_<dggs_type>``. Other parameters match
    ``h32cogp``.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    dggs2cogp(
        dggrid_spec(
            dggrid_instance, dggs_type, int(input_resolution), split_antimeridian, options
        ),
        input_path,
        output_path,
        agg_col=agg_col,
        id_col=dggrid_col or f"dggrid_{dggs_type.lower()}",
        min_res=min_res,
        source_resolution=int(input_resolution),
        verbose=verbose,
    )


def dggrid2cogp_cli():
    parser = cogp_arg_parser(
        "DGGRID", "dggrid_col", None, DEFAULT_MIN_RES, "dggrid_isea7h_10"
    )
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGRID_TYPES.keys(),
        help="DGGRID type",
    )
    parser.add_argument(
        "-ir",
        "--input_resolution",
        type=int,
        required=True,
        help="Resolution of the input DGGRID cells (SEQNUMs do not encode it).",
    )
    add_split_antimeridian_argument(parser)
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_from_cellids.",
    )
    args = parser.parse_args()
    run_cli(
        lambda: dggrid2cogp(
            create_dggrid_instance(),
            args.dggs_type,
            args.input,
            default_output_path(args),
            args.input_resolution,
            args.agg_col,
            args.id_col,
            args.min_res,
            split_antimeridian=args.split_antimeridian,
            options=json.loads(args.options) if args.options else None,
            verbose=args.verbose,
        )
    )


if __name__ == "__main__":
    dggrid2cogp_cli()

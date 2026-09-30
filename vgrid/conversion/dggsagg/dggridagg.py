"""
DGGRID Aggregate Module

Roll DGGRID cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

DGGRID SEQNUM IDs do not carry their resolution, so the input cells must all
be at ``input_resolution``. Each cell is assigned to the cell at
``resolution`` that contains its centroid, using two batched DGGRID calls.

Key Functions:
    dggrid_agg: Group DGGRID IDs by parent cell at a resolution
    dggridagg: Read cells into a GeoDataFrame and aggregate by parent
    dggridagg_cli: Command-line interface
"""

import json

import geopandas as gpd
from vgrid.utils.geometry import dggrid_num_edges, dggs_cell_row
from vgrid.utils.io import (
    create_dggrid_instance,
    prepare_compact_bags,
    validate_dggrid_resolution,
    validate_dggrid_type,
)
from vgrid.utils.constants import DGGRID_TYPES
from vgrid.conversion.dggs2geo.dggrid2geo import dggrid2geo
from vgrid.conversion.dggsagg.common import (
    agg_arg_parser,
    build_agg_rows,
    finish_agg_output,
    group_by_parent,
    print_structured,
)


def dggrid_parent_map(dggrid_instance, dggs_type, cell_ids, input_resolution, resolution):
    """Map each cell ID at ``input_resolution`` to its parent ID at ``resolution``."""
    cell_ids = [str(c) for c in cell_ids]
    if input_resolution == resolution:
        return {c: c for c in cell_ids}
    centroids = dggrid_instance.grid_cell_centroids_from_cellids(
        cell_ids, dggs_type, input_resolution
    )
    centroids = centroids[["global_id", "geometry"]].reset_index(drop=True)
    points = gpd.GeoDataFrame(
        {"geometry": centroids.geometry.values}, crs="EPSG:4326"
    )
    parents = dggrid_instance.cells_for_geo_points(
        geodf_points_wgs84=points,
        cell_ids_only=True,
        dggs_type=dggs_type,
        resolution=resolution,
    )
    parent_col = "name" if "name" in parents.columns else "seqnums"
    return {
        str(child): str(parent)
        for child, parent in zip(centroids["global_id"], parents[parent_col])
    }


def dggrid_agg(
    dggrid_instance,
    dggs_type,
    cell_ids,
    input_resolution,
    resolution,
    bags=None,
    verbose=True,
):
    """Group DGGRID cell IDs at ``input_resolution`` by their parent at ``resolution``."""
    dggs_type = validate_dggrid_type(dggs_type)
    input_resolution = validate_dggrid_resolution(dggs_type, int(input_resolution))
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))
    if resolution > input_resolution:
        raise ValueError(
            f"Parent resolution {resolution} is finer than input resolution "
            f"{input_resolution}."
        )
    parent_map = dggrid_parent_map(
        dggrid_instance, dggs_type, cell_ids, input_resolution, resolution
    )
    return group_by_parent(
        [str(c) for c in cell_ids],
        parent_map.__getitem__,
        "DGGRID",
        bags=bags,
        verbose=verbose,
    )


def dggridagg(
    dggrid_instance,
    dggs_type,
    input_data,
    input_resolution,
    resolution,
    dggrid_id=None,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    split_antimeridian=False,
    options=None,
    cell_metrics=False,
    verbose=True,
):
    """Aggregate DGGRID cells at ``input_resolution`` by parent at ``resolution``.

    The ID column defaults to ``dggrid_<dggs_type>``. Other parameters match
    ``h3agg``; with ``cell_metrics`` off the result is a DataFrame of the
    DGGRID id and the aggregated column.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    input_resolution = validate_dggrid_resolution(dggs_type, int(input_resolution))
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))
    id_col = f"dggrid_{dggs_type.lower()}"
    if dggrid_id is None:
        dggrid_id = id_col

    bags, agg_col = prepare_compact_bags(
        input_data,
        dggrid_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="DGGRID cells",
    )
    if bags is None:
        print(f"No DGGRID IDs found in <{dggrid_id}> field.")
        return None
    bags = {str(k): v for k, v in bags.items()}

    parent_ids = dggrid_agg(
        dggrid_instance,
        dggs_type,
        list(bags.keys()),
        input_resolution,
        resolution,
        bags=bags,
        verbose=verbose,
    )
    if not parent_ids:
        return None

    geoms = {}
    if cell_metrics:
        cells_gdf = dggrid2geo(
            dggrid_instance,
            dggs_type,
            parent_ids,
            resolution,
            split_antimeridian=split_antimeridian,
            options=options,
        )
        geoms = {
            str(cid): geom
            for cid, geom in zip(cells_gdf[id_col], cells_gdf.geometry)
        }
    num_edges = dggrid_num_edges(dggs_type)

    def row_fn(parent_id):
        return dggs_cell_row(
            id_col, parent_id, resolution, geoms[parent_id], num_edges, True
        )

    rows = build_agg_rows(
        parent_ids,
        bags,
        agg,
        agg_col,
        dggrid_id,
        row_fn,
        cell_metrics,
        "DGGRID",
        verbose=verbose,
    )
    return finish_agg_output(
        rows, cell_metrics, output_format, input_data, f"dggrid_{dggs_type.lower()}"
    )


def dggridagg_cli():
    """Command-line interface for dggridagg."""
    parser = agg_arg_parser("DGGRID", "dggridcompact")
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
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable Antimeridian splitting",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_from_cellids.",
    )
    args = parser.parse_args()
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error parsing options JSON: {e}")
            return
    result = dggridagg(
        create_dggrid_instance(),
        args.dggs_type,
        args.input,
        input_resolution=args.input_resolution,
        resolution=args.resolution,
        dggrid_id=args.cellid,
        agg=args.agg,
        numeric_col=args.numeric_col,
        output_format=args.output_format,
        split_antimeridian=args.split_antimeridian,
        options=options,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    print_structured(args.output_format, result)

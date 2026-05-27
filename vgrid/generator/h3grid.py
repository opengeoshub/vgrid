# Reference: https://observablehq.com/@claude-ducharme/h3-map
# https://h3-snow.streamlit.app/

import argparse
from shapely.geometry import box, shape
from shapely.ops import unary_union
from tqdm import tqdm
import geopandas as gpd
import h3
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS, MAX_CELLS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    is_full_world_bbox,
    validate_bbox,
    validate_h3_resolution,
    convert_to_output_format,
)
from vgrid.conversion.dggs2geo.h32geo import h32geo


def _polygon_to_h3_cells_exprimental(polygon, resolution, contain=None):
    """
    Return H3 cell IDs covering a Shapely polygon.

    Uses ``geo_to_cells`` by default. When ``contain`` is set, uses
    ``polygon_to_cells_experimental``, which requires H3 ``LatLngPoly`` input.
    """
    if contain is None:
        return h3.geo_to_cells(polygon, resolution)
    outer = [(lat, lon) for lon, lat in polygon.exterior.coords[:-1]]
    holes = [
        [(lat, lon) for lon, lat in ring.coords[:-1]]
        for ring in polygon.interiors
    ]
    h3_poly = h3.LatLngPoly(outer, *holes)
    return h3.polygon_to_cells_experimental(h3_poly, resolution, contain=contain)


def h3_grid(resolution, fix_antimeridian=None):
    resolution = validate_h3_resolution(resolution)
    total_cells = h3.get_num_cells(resolution)
    if total_cells > MAX_CELLS:
        raise ValueError(
            f"Resolution {resolution} will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
        )
    else:
        base_cells = h3.get_res0_cells()
        h3_records = []
        # Progress bar for base cells
        with tqdm(total=total_cells, desc="Generating H3 DGGS", unit=" cells") as pbar:
            for cell in base_cells:
                child_cells = h3.cell_to_children(cell, resolution)
                # Progress bar for child cells
                for child_cell in child_cells:
                    cell_polygon = h32geo(child_cell, fix_antimeridian=fix_antimeridian)
                    h3_id = str(child_cell)
                    num_edges = 6
                    if h3.is_pentagon(h3_id):
                        num_edges = 5
                    record = geodesic_dggs_to_geoseries(
                        "h3", h3_id, resolution, cell_polygon, num_edges
                    )
                    h3_records.append(record)
                    pbar.update(1)

        return gpd.GeoDataFrame(h3_records, geometry="geometry", crs="EPSG:4326")


def h3_grid_within_bbox(resolution, bbox, fix_antimeridian=None):
    resolution = validate_h3_resolution(resolution)
    bbox = validate_bbox(bbox)
    if is_full_world_bbox(bbox):
        return h3_grid(resolution, fix_antimeridian=fix_antimeridian)

    bbox_polygon = box(*bbox)
    bbox_cells = h3.geo_to_cells(bbox_polygon, resolution)
    # bbox_cells = _polygon_to_h3_cells_exprimental(bbox_polygon, resolution, contain="overlap")
    total_cells = len(bbox_cells)
    if total_cells > MAX_CELLS:
        raise ValueError(
            f"Resolution {resolution} within bounding box {bbox} will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
        )
    else:
        h3_records = []
        for cell_id in tqdm(bbox_cells, desc="Generating H3 DGGS"):
            cell_polygon = h32geo(cell_id, fix_antimeridian=fix_antimeridian)
            if cell_polygon.intersects(bbox_polygon):
                h3_id = str(cell_id)
                num_edges = 6
                if h3.is_pentagon(h3_id):
                    num_edges = 5
                record = geodesic_dggs_to_geoseries(
                    "h3", h3_id, resolution, cell_polygon, num_edges
                )
                h3_records.append(record)

        return gpd.GeoDataFrame(h3_records, geometry="geometry", crs="EPSG:4326")


def h3_grid_ids(resolution, fix_antimeridian=None):
    """
    Generate a list of H3 cell IDs for the whole world at a given resolution.

    Args:
        resolution (int): H3 resolution [0..15]

    Returns:
        list[str]: List of H3 cell IDs as strings
    """
    resolution = validate_h3_resolution(resolution)
    total_cells = h3.get_num_cells(resolution)
    base_cells = h3.get_res0_cells()
    h3_ids = []
    with tqdm(total=total_cells, desc="Generating H3 IDs", unit=" cells") as pbar:
        for cell in base_cells:
            child_cells = h3.cell_to_children(cell, resolution)
            for child_cell in child_cells:
                h3_ids.append(str(child_cell))
                pbar.update(1)

    return h3_ids


def h3_grid_within_bbox_ids(resolution, bbox, fix_antimeridian=None):
    """
    Generate a list of H3 cell IDs that intersect a bounding box.

    Args:
        resolution (int): H3 resolution [0..15]
        bbox (list[float]): [min_lon, min_lat, max_lon, max_lat]

    Returns:
        list[str]: List of H3 cell IDs as strings that intersect the bbox
    """
    resolution = validate_h3_resolution(resolution)
    bbox = validate_bbox(bbox)
    if is_full_world_bbox(bbox):
        return h3_grid_ids(resolution, fix_antimeridian=fix_antimeridian)

    bbox_polygon = box(*bbox)
    bbox_cells = h3.geo_to_cells(bbox_polygon, resolution)
    total_cells = len(bbox_cells)
    h3_ids = []
    for cell_id in tqdm(bbox_cells, total=total_cells, desc="Generating H3 IDs"):
        cell_polygon = h32geo(cell_id, fix_antimeridian=fix_antimeridian)
        if cell_polygon.intersects(bbox_polygon):
            h3_ids.append(str(cell_id))

    return h3_ids


def h3grid(resolution, bbox=None, output_format="gpd", fix_antimeridian=None):
    """
    Generate H3 grid for pure Python usage.

    Args:
        resolution (int): H3 resolution [0..15]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format handled entirely by convert_to_output_format

    Returns:
        Delegated to convert_to_output_format
    """
    if bbox is None:
        h3_gdf = h3_grid(resolution, fix_antimeridian=fix_antimeridian)
    else:
        h3_gdf = h3_grid_within_bbox(
            resolution, bbox, fix_antimeridian=fix_antimeridian
        )
    output_name = f"h3_grid_{resolution}"
    return convert_to_output_format(h3_gdf, output_format, output_name)


def h3grid_cli():
    """CLI interface for generating H3 grid."""
    parser = argparse.ArgumentParser(description="Generate H3 DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..15]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box: <min_lon> <min_lat> <max_lon> <max_lat> (default is the whole world)",
    )
    parser.add_argument(
        "-f", "--output_format", type=str, choices=OUTPUT_FORMATS, default="gpd"
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )

    args = parser.parse_args()
    try:
        result = h3grid(
            args.resolution,
            args.bbox,
            args.output_format,
            fix_antimeridian=args.fix_antimeridian,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return

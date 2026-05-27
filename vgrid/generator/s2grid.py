"""
S2 Grid Generator Module

Generates S2 DGGS grids for specified resolutions and bounding boxes with automatic cell generation and validation.

Key Functions:
- s2_grid(): Main grid generation function with bounding box support
- s2_grid_ids(): Returns list of S2 cell tokens for given resolution and bbox
- s2grid(): User-facing function with multiple output formats
- s2grid_cli(): Command-line interface for grid generation

Reference:
    https://github.com/aaliddell/s2cell,
    https://medium.com/@claude.ducharme/selecting-a-geo-representation-81afeaf3bf01
    https://github.com/sidewalklabs/s2
    https://github.com/google/s2geometry/tree/master/src/python
    https://github.com/google/s2geometry
    https://gis.stackexchange.com/questions/293716/creating-shapefile-of-s2-cells-for-given-level
    https://s2.readthedocs.io/en/latest/quickstart.html

"""

import argparse
from shapely.geometry import shape
from shapely.ops import unary_union
import geopandas as gpd
from tqdm import tqdm
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.dggs import s2
from vgrid.utils.io import validate_bbox, validate_s2_resolution, convert_to_output_format
from vgrid.conversion.dggs2geo.s22geo import s22geo


def s2_grid(resolution, bbox, fix_antimeridian=None):
    """
    Generate an S2 DGGS grid for a given resolution and bounding box.
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
    Args:
        resolution (int): S2 level [0..30]
        bbox (list[float]): [min_lon, min_lat, max_lon, max_lat]
        fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
    Returns:
        geopandas.GeoDataFrame: GeoDataFrame containing the S2 DGGS grid
    """
    resolution = validate_s2_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    level = resolution
    cell_ids = []
    coverer = s2.RegionCoverer()
    coverer.min_level = level
    coverer.max_level = level
    # coverer.max_cells = 1000_000  # Adjust as needed
    # coverer.max_cells = 0  # Adjust as needed

    region = s2.LatLngRect(
        s2.LatLng.from_degrees(min_lat, min_lon),
        s2.LatLng.from_degrees(max_lat, max_lon),
    )

    # Get the covering cells
    covering = coverer.get_covering(region)

    # Convert the covering cells to S2 cell IDs
    for cell_id in covering:
        cell_ids.append(cell_id)

    s2_rows = []
    num_edges = 4

    for cell_id in tqdm(cell_ids, desc="Generating DGGS", unit=" cells"):
        # Generate a Shapely Polygon
        cell_polygon = s22geo(cell_id.to_token(), fix_antimeridian=fix_antimeridian)
        s2_token = cell_id.to_token()
        row = geodesic_dggs_to_geoseries(
            "s2", s2_token, resolution, cell_polygon, num_edges
        )
        s2_rows.append(row)

    return gpd.GeoDataFrame(s2_rows, geometry="geometry", crs="EPSG:4326")


def s2_grid_ids(resolution, bbox, fix_antimeridian=None):
    """
    Return a list of S2 cell tokens for a given resolution and bounding box.
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
    Args:
        resolution (int): S2 level [0..30]
        bbox (list[float]): [min_lon, min_lat, max_lon, max_lat]
        fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
    Returns:
        list[str]: List of S2 cell tokens
    """
    resolution = validate_s2_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    level = resolution
    coverer = s2.RegionCoverer()
    coverer.min_level = level
    coverer.max_level = level
    region = s2.LatLngRect(
        s2.LatLng.from_degrees(min_lat, min_lon),
        s2.LatLng.from_degrees(max_lat, max_lon),
    )
    covering = coverer.get_covering(region)
    return [cell_id.to_token() for cell_id in covering]


def s2grid(resolution, bbox=None, output_format="gpd", fix_antimeridian=None):
    """
    Generate S2 grid for pure Python usage.

    Args:
        resolution (int): S2 resolution [0..30]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output output_format ('geojson', 'csv', etc.). Defaults to None (list of S2 tokens).
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none. Defaults to None.
    Returns:
        dict or list: GeoJSON FeatureCollection, list of S2 tokens, or file path depending on output_format
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        num_cells = 6 * (4**resolution)
        if num_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {num_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
    gdf = s2_grid(resolution, bbox, fix_antimeridian=fix_antimeridian)
    output_name = f"s2_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def s2grid_cli():
    """CLI interface for generating S2 grid."""
    parser = argparse.ArgumentParser(description="Generate S2 DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..30]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the output_format: min_lon min_lat max_lon max_lat (default is the whole world)",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
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
        result = s2grid(
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


if __name__ == "__main__":
    s2grid_cli()

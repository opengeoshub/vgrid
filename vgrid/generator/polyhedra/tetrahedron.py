"""
Tetrahedron Grid Generator Module

Generates Tetrahedron DGGS grids with automatic cell generation and validation for tetrahedral polyhedron representation.

Key Functions:
- tetrahedron(): Main grid generation function for tetrahedral polyhedron
- tetrahedron_cli(): Command-line interface for grid generation
"""

import argparse
from shapely.geometry import Polygon
import geopandas as gpd
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.io import convert_to_output_format


def constructGeometry(facet):
    """Create a Polygon with the vertices (longitude, latitude)."""
    poly = Polygon([(v[0], v[1]) for v in facet])  # (lon, lat)
    return poly


def tetrahedron(output_format="gpd"):
    """
    Generate tetrahedron faces as GeoDataFrame.

    Args:
        output_format (str, optional): Output format ('geojson', 'csv', etc.). Defaults to "gpd".

    Returns:
        dict or list: GeoJSON FeatureCollection, list of tokens, or file path depending on output_format
    """

    # Define facets with coordinates and cell_ids
    facets = [
        {
            "cell_id": "0",
            "coordinates": [
                [-180.0, 0.0],
                [-180.0, 90.0],
                [-90.0, 90.0],
                [0.0, 90.0],
                [0.0, 0.0],
                [-90.0, 0.0],
                [-180.0, 0.0],
            ],
        },
        {
            "cell_id": "1",
            "coordinates": [
                [0.0, 0.0],
                [0.0, 90.0],
                [90.0, 90.0],
                [180.0, 90.0],
                [180.0, 0.0],
                [90.0, 0.0],
                [0.0, 0.0],
            ],
        },
        {
            "cell_id": "2",
            "coordinates": [
                [-180.0, -90.0],
                [-180.0, 0.0],
                [-90.0, 0.0],
                [0.0, 0.0],
                [0.0, -90.0],
                [-90.0, -90.0],
                [-180.0, -90.0],
            ],
        },
        {
            "cell_id": "3",
            "coordinates": [
                [0.0, -90.0],
                [0.0, 0.0],
                [90.0, 0.0],
                [180.0, 0.0],
                [180.0, -90.0],
                [90.0, -90.0],
                [0.0, -90.0],
            ],
        },
    ]

    # Create lists to store geometries and properties
    geometries = []
    zone_ids = []

    for facet in facets:
        geometry = Polygon(facet["coordinates"])
        geometries.append(geometry)
        zone_ids.append(facet["cell_id"])

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame({"cell_id": zone_ids}, geometry=geometries, crs="EPSG:4326")

    output_name = "tetrahedron"
    return convert_to_output_format(gdf, output_format, output_name)


def tetrahedron_cli():
    """CLI interface for generating tetrahedron faces."""
    parser = argparse.ArgumentParser(description="Generate tetrahedron faces.")

    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    args = parser.parse_args()
    try:
        result = tetrahedron(output_format=args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    tetrahedron_cli()

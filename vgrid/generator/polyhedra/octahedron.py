"""
Octahedron Grid Generator Module

Generates Octahedron DGGS grids with automatic cell generation and validation for octahedral polyhedron representation.

Key Functions:
- octahedron(): Main grid generation function for octahedral polyhedron
- octahedron_cli(): Command-line interface for grid generation
"""

import argparse
from shapely.geometry import Polygon, LinearRing
import geopandas as gpd
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.io import convert_to_output_format


def constructGeometry(facet):
    """Construct geometry from facet coordinates."""
    vertexTuples = facet[:4]
    # Create a LinearRing with the vertices
    ring = LinearRing(
        [(vT[1], vT[0]) for vT in vertexTuples]
    )  # sequence: lon, lat (x,y)

    # Create a Polygon from the LinearRing
    poly = Polygon(ring)
    return poly


def octahedron(output_format="gpd"):
    """
    Generate octahedron faces as GeoDataFrame.

    Args:
        output_format (str, optional): Output format ('geojson', 'csv', etc.). Defaults to "gpd".

    Returns:
        dict or list: GeoJSON FeatureCollection, list of tokens, or file path depending on output_format
    """
    # Define coordinate points for octahedron faces
    p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (
        (90.0, -180.0),
        (90.0, -90.0),
        (90.0, 0.0),
        (90.0, 90.0),
        (90.0, 180.0),
    )
    p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (
        (0.0, -180.0),
        (0.0, -90.0),
        (0.0, 0.0),
        (0.0, 90.0),
        (0.0, 180.0),
    )
    n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (
        (-90.0, -180.0),
        (-90.0, -90.0),
        (-90.0, 0.0),
        (-90.0, 90.0),
        (-90.0, 180.0),
    )

    # Define initial facets for octahedron
    initial_facets = [
        [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
        [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
        [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
        [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
        [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
        [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
        [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
        [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
    ]

    # Create lists to store geometries and properties
    geometries = []
    cell_ids = []

    for i, facet in enumerate(initial_facets):
        geometry = constructGeometry(facet)
        geometries.append(geometry)
        cell_ids.append(str(i + 1))

    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame({"cell_id": cell_ids}, geometry=geometries, crs="EPSG:4326")

    output_name = "octahedron"
    return convert_to_output_format(gdf, output_format, output_name)


def octahedron_cli():
    """CLI interface for generating octahedron faces."""
    parser = argparse.ArgumentParser(description="Generate octahedron faces.")

    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    args = parser.parse_args()
    try:
        result = octahedron(output_format=args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    octahedron_cli()

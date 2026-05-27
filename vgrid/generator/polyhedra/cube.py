"""
Cube S2 Grid Generator Module.
This module provides functionality to generate S2 cube.
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
from vgrid.generator.s2grid import s2_grid
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.io import convert_to_output_format


def cube(output_format="gpd"):
    """
    Generate S2 grid for pure Python usage.

    Args:
        resolution (int): S2 resolution [0..30]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output output_format ('geojson', 'csv', etc.). Defaults to None (list of S2 tokens).

    Returns:
        dict or list: GeoJSON FeatureCollection, list of S2 tokens, or file path depending on output_format
    """
    bbox = [-180, -90, 180, 90]
    gdf = s2_grid(0, bbox)
    output_name = "cube"
    return convert_to_output_format(gdf, output_format, output_name)


def cube_cli():
    """CLI interface for generating S2 cube."""
    parser = argparse.ArgumentParser(description="Generate S2 DGGS.")

    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    args = parser.parse_args()
    try:
        result = cube(resolution=0, output_format=args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    cube_cli()

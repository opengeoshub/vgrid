"""
Dodecahedron Grid Generator Module

Generates Dodecahedron DGGS grids with automatic cell generation and validation for dodecahedral polyhedron representation using DGGAL library.

Key Functions:
- dodecahedron(): Main grid generation function for dodecahedral polyhedron
- dodecahedron_cli(): Command-line interface for grid generation
"""

import sys
import argparse
import geopandas as gpd
from dggal import *
from tqdm import tqdm
from vgrid.utils.constants import DGGAL_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.io import dggal_convert_to_output_format
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo

app = Application(appGlobals=globals())
pydggal_setup(app)


def dodecahedron(output_format="gpd", split_antimeridian=False):
    """
    Generate a DGGAL grid using the dggal library directly.

    When output_format is provided, save to the current folder using a predefined
    name (e.g., "<dggs_type>_grid_<resolution>.*"), mirroring h3grid behavior.
    Returns either a GeoDataFrame, a path/string, a dict, or a list depending on output_format.

    Args:
        output_format (str, optional): Output format.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
    """
    dggs_type = "isea3h"
    resolution = 0
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()
    # Set up bbox for listZones if provided
    geo_extent = wholeWorld  # Default to whole world
    # Call listZones to get all zones at the specified resolution
    zones = dggrs.listZones(resolution, geo_extent)

    dggal_records = []

    for zone in tqdm(zones, desc="Generating dodecahedron"):
        try:
            zone_id = dggrs.getZoneTextID(zone)
            zone_resolution = dggrs.getZoneLevel(zone)
            num_edges = dggrs.countZoneEdges(zone)

            # Convert zone to geometry using dggal2geo
            cell_polygon = dggal2geo(
                dggs_type, zone_id, split_antimeridian=split_antimeridian
            )

            # Create record using geodesic_dggs_to_geoseries
            record = geodesic_dggs_to_geoseries(
                f"dggal_{dggs_type}", zone_id, zone_resolution, cell_polygon, num_edges
            )
            dggal_records.append(record)

        except Exception as e:
            print(f"Error processing zone: {e}", file=sys.stderr)
            continue

    # Create GeoDataFrame from records
    gdf = gpd.GeoDataFrame(dggal_records, geometry="geometry", crs="EPSG:4326")
    base_name = "dodecahedron"

    return dggal_convert_to_output_format(gdf, output_format, output_name=base_name)


def dodecahedron_cli():
    """CLI interface for generating dodecahedron."""
    parser = argparse.ArgumentParser(description="Generate dodecahedron.")

    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    args = parser.parse_args()
    try:
        result = dodecahedron(output_format=args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    dodecahedron_cli()

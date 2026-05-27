"""
DGGAL Grid Generator Module

Generates DGGAL (Discrete Global Grids with Adaptive Localization) grids for multiple grid types with automatic cell generation and validation.

Key Functions:
- dggalgen(): User-facing function with multiple output formats
- dggalgen_cli(): Command-line interface for grid generation
"""

import argparse
import sys
import geopandas as gpd
from vgrid.utils.io import dggal_convert_to_output_format
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS, DGGAL_TYPES
from vgrid.utils.io import validate_bbox, validate_dggal_resolution, validate_dggal_type
from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo

from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from tqdm import tqdm

# Import dggal library
from dggal import *

# Initialize dggal application
app = Application(appGlobals=globals())
pydggal_setup(app)


def dggalgen(
    dggs_type: str = "gnosis",
    resolution: int = 1,
    bbox: tuple[float, float, float, float] | None = None,
    compact: bool = False,
    output_format: str | None = None,
    split_antimeridian: bool = False,
):
    """
    Generate a DGGAL grid using the dggal library directly.

    When output_format is provided, save to the current folder using a predefined
    name (e.g., "<dggs_type>_grid_<resolution>.*"), mirroring h3grid behavior.
    Returns either a GeoDataFrame, a path/string, a dict, or a list depending on output_format.

    Parameters
    ----------
    dggs_type : str, default "gnosis"
        DGGAL DGGS type.
    resolution : int, default 1
        Resolution level.
    bbox : tuple[float, float, float, float] | None, optional
        Bounding box as (min_lon, min_lat, max_lon, max_lat).
    compact : bool, default False
        Whether to compact zones.
    output_format : str | None, optional
        Output format.
    split_antimeridian : bool, default False
        When True, apply antimeridian fixing to the resulting polygons.
    """

    # Validate resolution against per-type bounds
    dggs_type = validate_dggal_type(dggs_type)
    resolution = validate_dggal_resolution(dggs_type, resolution)

    # Create the appropriate DGGS instance
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()

    # Set up bbox for listZones if provided
    geo_extent = wholeWorld  # Default to whole world
    if bbox:
        # bbox should be (min_lon, min_lat, max_lon, max_lat)
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
        # Validate bbox coordinates
        valid_lat = (-90 < min_lat < 90) and (-90 < max_lat < 90)
        if valid_lat:
            ll = GeoPoint(min_lat, min_lon)
            ur = GeoPoint(max_lat, max_lon)
            geo_extent = GeoExtent(ll, ur)
        else:
            print(
                "Invalid bounding box coordinates, using whole world", file=sys.stderr
            )
            geo_extent = wholeWorld
    # Call listZones to get all zones at the specified resolution
    zones = dggrs.listZones(resolution, geo_extent)
    # Check if zones is None or empty before compacting
    if zones is None or len(zones) == 0:
        print("No zones found", file=sys.stderr)
        return None

    # Only compact if requested and zones exist
    if compact:
        compacted_zones = dggrs.compactZones(zones)
        if compacted_zones is not None:
            zones = compacted_zones

    # Zones found successfully, proceed with processing

    # Process zones directly with tqdm progress bar
    dggal_records = []
    options = {}

    for zone in tqdm(zones, desc=f"Generating {dggs_type.upper()} DGGS"):
        try:
            zone_id = dggrs.getZoneTextID(zone)
            zone_resolution = dggrs.getZoneLevel(zone)
            num_edges = dggrs.countZoneEdges(zone)

            # Convert zone to geometry using dggal2geo
            cell_polygon = dggal2geo(
                dggs_type, zone_id, options, split_antimeridian=split_antimeridian
            )

            # Create record using geodesic_dggs_to_geoseries
            record = geodesic_dggs_to_geoseries(
                f"dggal_{dggs_type}", zone_id, zone_resolution, cell_polygon, num_edges
            )
            dggal_records.append(record)

        except Exception as e:
            print(f"Error processing zone: {e}", file=sys.stderr)
            continue

    if dggal_records:
        # Create GeoDataFrame from records
        gdf = gpd.GeoDataFrame(dggal_records, geometry="geometry", crs="EPSG:4326")
        base_name = f"{dggs_type}_grid_{resolution}"
        return dggal_convert_to_output_format(
            gdf, output_format, output_name=base_name
        )
    else:
        print("No valid zones found for the specified parameters.", file=sys.stderr)
        return None


def dggalgen_cli():
    parser = argparse.ArgumentParser(description="Generate grid via dggal library.")
    parser.add_argument(
        "-dggs", "--dggs_type", required=True, type=str, choices=DGGAL_TYPES.keys()
    )
    parser.add_argument(
        "-r",
        "--resolution",
        required=True,
        type=int,
        help="Resolution (integer)",
    )
    parser.add_argument(
        "-c",
        "--compact",
        dest="compact",
        required=False,
        action="store_true",
        help="Compact zones",
    )
    parser.add_argument(
        "-bbox",
        "--bbox",
        dest="bbox",
        required=False,
        type=str,
        help="Bounding box as 'min_lon,min_lat,max_lon,max_lat' (e.g., '100,10,110,20')",
    )
    parser.add_argument(
        "-f", "--output_format", type=str, default="gpd", choices=OUTPUT_FORMATS
    )
    # No custom output path; files are saved in current folder with predefined names
    args = parser.parse_args()

    # Parse bbox if provided
    bbox_tuple = None
    if args.bbox:
        try:
            bbox_parts = args.bbox.split(",")
            if len(bbox_parts) == 4:
                bbox_tuple = tuple(float(x) for x in bbox_parts)
            else:
                print(
                    "Error: bbox must be in format 'min_lon,min_lat,max_lon,max_lat'",
                    file=sys.stderr,
                )
                sys.exit(1)
        except ValueError:
            print("Error: bbox coordinates must be numeric", file=sys.stderr)
            sys.exit(1)

    result = dggalgen(
        dggs_type=args.dggs_type,
        resolution=args.resolution,
        bbox=bbox_tuple,
        compact=args.compact,
    )
    if result is None:
        sys.exit(1)
    # Structured formats print to console; file outputs are announced by the IO helper
    if args.output_format in STRUCTURED_FORMATS:
        print(result)
    sys.exit(0)


if __name__ == "__main__":
    dggalgen_cli()

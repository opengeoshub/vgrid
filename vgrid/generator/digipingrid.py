"""
DIGIPIN Grid Generator Module

Generates DIGIPIN DGGS grids for specified resolutions with automatic cell generation and validation for India region.

Key Functions:
- digipin_grid(): Main grid generation function for India region
- digipin_grid_within_bbox(): Grid generation within bounding box
- digipingrid(): User-facing function with multiple output formats
- digipingrid_cli(): Command-line interface for grid generation

Note: DIGIPIN is a geocoding system for India with bounds (lat: 2.5-38.5, lon: 63.5-99.5)
"""

import argparse
from shapely.geometry import Polygon, shape
from shapely.ops import unary_union
from tqdm import tqdm
from vgrid.utils.constants import (
    MAX_CELLS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)
from vgrid.utils.geometry import graticule_dggs_to_geoseries
import geopandas as gpd
from vgrid.dggs.digipin import BOUNDS
from vgrid.conversion.latlon2dggs import latlon2digipin
from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo
from vgrid.utils.io import (
    validate_bbox,
    validate_digipin_resolution,
    convert_to_output_format,
)
from vgrid.conversion.dggscompact.digipincompact import digipin_compact
from vgrid.dggs.digipin import digipin_resolution


def _digipin_row_from_id(digipin_code):
    cell_polygon = digipin2geo(digipin_code)
    if isinstance(cell_polygon, str):
        raise ValueError(f"Invalid DIGIPIN cell: {digipin_code}")
    cell_resolution = digipin_resolution(digipin_code)
    return graticule_dggs_to_geoseries(
        "digipin", digipin_code, cell_resolution, cell_polygon
    )


def digipin_grid(resolution, bbox=None, compact=False):
    """
    Generate DIGIPIN grid at the given resolution.

    Parameters
    ----------
    resolution : int
        DIGIPIN resolution level [1-10]
    bbox : list, optional
        Bounding box [min_lon, min_lat, max_lon, max_lat].
        If None, defaults to entire India region.
    compact : bool, optional
        Enable DIGIPIN compact mode to reduce cell count.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame containing DIGIPIN cells with geometries and metadata
    """
    digipin_ids = digipin_grid_ids(resolution, bbox=bbox)
    if compact:
        digipin_ids = digipin_compact(digipin_ids)

    digipin_records = []
    for digipin_code in tqdm(
        digipin_ids, desc="Generating DIGIPIN DGGS", unit=" cells"
    ):
        try:
            digipin_records.append(_digipin_row_from_id(digipin_code))
        except Exception:
            continue

    return gpd.GeoDataFrame(digipin_records, geometry="geometry", crs="EPSG:4326")


def digipin_grid_ids(resolution, bbox=None, compact=False):
    """
    Return a list of DIGIPIN IDs at the given resolution.

    Parameters
    ----------
    resolution : int
        DIGIPIN resolution level [1-10]
    bbox : list, optional
        Bounding box [min_lon, min_lat, max_lon, max_lat].
        If None, defaults to entire India region.
    compact : bool, optional
        Enable DIGIPIN compact mode to reduce cell count.

    Returns
    -------
    list
        List of DIGIPIN cell IDs
    """
    resolution = validate_digipin_resolution(resolution)

    # Default to India bounds if no bbox provided
    if bbox is None:
        bbox = [BOUNDS["minLon"], BOUNDS["minLat"], BOUNDS["maxLon"], BOUNDS["maxLat"]]

    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)

    # Constrain to DIGIPIN bounds (India region)
    min_lat = max(min_lat, BOUNDS["minLat"])
    min_lon = max(min_lon, BOUNDS["minLon"])
    max_lat = min(max_lat, BOUNDS["maxLat"])
    max_lon = min(max_lon, BOUNDS["maxLon"])

    # Calculate sampling density based on resolution
    base_width = 9.0  # degrees at resolution 1
    factor = 0.25 ** (resolution - 1)  # each level divides by 4
    sample_width = base_width * factor

    seen_cells = set()
    ids = []

    # Sample points across the bounding box
    lon = min_lon
    while lon <= max_lon:
        lat = min_lat
        while lat <= max_lat:
            try:
                # Get DIGIPIN code for this point at the specified resolution
                digipin_code = latlon2digipin(lat, lon, resolution)

                if digipin_code == "Out of Bound":
                    lat += sample_width
                    continue

                if digipin_code not in seen_cells:
                    seen_cells.add(digipin_code)
                    ids.append(digipin_code)

            except Exception:
                # Skip cells with errors
                pass

            lat += sample_width
        lon += sample_width

    if compact:
        ids = digipin_compact(ids)
    return ids


def digipingrid(resolution, bbox=None, output_format="gpd", compact=False):
    """
    Generate DIGIPIN grid for pure Python usage.

    Parameters
    ----------
    resolution : int
        DIGIPIN resolution level [1-10]
    bbox : list, optional
        Bounding box [min_lon, min_lat, max_lon, max_lat].
        Defaults to None (entire India region)
    output_format : str, optional
        Output format ('geojson', 'csv', 'gpd', 'shapefile', 'gpkg', 'parquet',
        or None for list of DIGIPIN IDs). Defaults to 'gpd'.
    compact : bool, optional
        Enable DIGIPIN compact mode to reduce cell count.

    Returns
    -------
    dict, list, or str
        Output in the requested format or file path

    Examples
    --------
    >>> # Generate grid for entire India at resolution 3
    >>> gdf = digipingrid(3)

    >>> # Generate grid for a specific region
    >>> bbox = [77.0, 28.0, 78.0, 29.0]  # Delhi region
    >>> gdf = digipingrid(5, bbox=bbox)

    >>> # Generate and save as GeoJSON
    >>> result = digipingrid(4, output_format='geojson')
    """
    # Rough estimate: 4^resolution cells for the entire region
    if bbox is None:
        total_cells = (
            4**resolution * 4
        )  # Approximate, as India is ~36°x36° = 4 base cells
        if total_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate approximately {total_cells} cells "
                f"which exceeds the limit of {MAX_CELLS}"
            )

    gdf = digipin_grid(resolution, bbox=bbox, compact=compact)

    output_name = f"digipin_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def digipingrid_cli():
    """
    Command-line interface for DIGIPIN grid generation.
    """
    parser = argparse.ArgumentParser(
        description="Generate DIGIPIN DGGS grid for India region."
    )
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [1..10]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is India region)",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable DIGIPIN compact mode to reduce cell count",
    )
    args = parser.parse_args()

    try:
        result = digipingrid(
            args.resolution, args.bbox, args.output_format, compact=args.compact
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    digipingrid_cli()

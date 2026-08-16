"""
VGRID Grid Generator Module

Generates VGRID DGGS grids for specified resolutions and bounding boxes with automatic cell generation and validation.

Key Functions:
- vgrid_gen_ids(): Generate VGRID IDs only (returns list of strings)
- vgrid_gen(): Main grid generation function with bounding box support (returns GeoDataFrame)
- vgridgen(): User-facing function with multiple output formats
- vgridgen_cli(): Command-line interface for grid generation
"""

import argparse
import geopandas as gpd
from shapely.geometry import Polygon
from tqdm import tqdm
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import (
    convert_to_output_format,
    validate_bbox,
    validate_vgrid_resolution,
)
from vgrid.dggs.vgrid import VGRID


def vgrid_gen_ids(
    vgrid_instance: VGRID, resolution: int = 0, bbox: list[float] = None
) -> list[str]:
    """
    Generate VGRID IDs for a given VGRID instance and bounding box.

    Args:
        vgrid_instance (VGRID): VGRID instance
        resolution (int, optional): Resolution level to generate. If None, uses VGRID instance resolution.
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).

    Returns:
        list[str]: List of VGRID cell IDs within the bounding box
    """
    if bbox is None:
        # Default to whole world
        bbox = [-180, -90, 180, 90]
    else:
        bbox = validate_bbox(bbox)

    # Use specified resolution or VGRID instance resolution
    resolution = validate_vgrid_resolution(resolution)

    # Calculate cell size for the target resolution
    # Use the same logic as tiles_for_bounding_box in VGRID.py
    # Generate tiles for the bounding box
    tiles = vgrid_instance.tiles_for_bounding_box(
        bbox[0], bbox[1], bbox[2], bbox[3], resolution
    )

    vgrid_ids = []

    with tqdm(
        total=len(tiles),
        desc=f"Generating VGRID IDs resolution {resolution}",
        unit=" cells",
    ) as pbar:
        for res, tile_index in tiles:
            # Get the VGRID ID for this tile
            vgrid_id = vgrid_instance.to_vgrid_id(resolution, tile_index, 0)
            vgrid_ids.append(str(vgrid_id))
            pbar.update(1)

    return vgrid_ids


def vgrid_gen(
    vgrid_instance: VGRID, resolution: int = 0, bbox: list[float] = None
) -> gpd.GeoDataFrame:
    """
    Generate a VGRID grid for a given VGRID instance and bounding box.

    Args:
        vgrid_instance (VGRID): VGRID instance
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        resolution (int, optional): Resolution level to generate. If None, uses VGRID instance resolution.

    Returns:
        GeoDataFrame: VGRID grid cells within the bounding box
    """
    if bbox is None:
        # Default to whole world
        bbox = [-180, -90, 180, 90]
    else:
        bbox = validate_bbox(bbox)

    # Use specified resolution or VGRID instance resolution
    resolution = validate_vgrid_resolution(resolution)

    # Calculate cell size for the target resolution
    # Use the same logic as tiles_for_bounding_box in VGRID.py
    base_cell_size = vgrid_instance.cell_size
    if resolution >= 0:
        # For positive resolutions, divide both longitude and latitude cell sizes by sqrt(aperture)^resolution
        # This ensures total cells increase by aperture factor, not aperture^2
        factor = (vgrid_instance.aperture**0.5) ** resolution
        target_cell_size_lon = base_cell_size / factor
        target_cell_size_lat = base_cell_size / factor
    else:
        # For negative resolutions, multiply both longitude and latitude cell sizes by sqrt(aperture)^|resolution|
        factor = (vgrid_instance.aperture**0.5) ** abs(resolution)
        target_cell_size_lon = base_cell_size * factor
        target_cell_size_lat = base_cell_size * factor
    target_total_columns = int(360 / target_cell_size_lon)

    # Generate tiles for the bounding box
    tiles = vgrid_instance.tiles_for_bounding_box(
        bbox[0], bbox[1], bbox[2], bbox[3], resolution
    )

    vgrid_records = []

    with tqdm(
        total=len(tiles),
        desc=f"Generating VGRID resolution {resolution}",
        unit=" cells",
    ) as pbar:
        for res, tile_index in tiles:
            # Get the bottom-left corner of the tile
            vgrid_id = vgrid_instance.to_vgrid_id(resolution, tile_index, 0)

            # Calculate lat/lon using target resolution's cell sizes
            lat = (tile_index // target_total_columns) * target_cell_size_lat - 90
            lon = (tile_index % target_total_columns) * target_cell_size_lon - 180

            # Create the tile polygon with proper bounds
            cell_polygon = Polygon(
                [
                    (lon, lat),
                    (min(180, lon + target_cell_size_lon), lat),
                    (
                        min(180, lon + target_cell_size_lon),
                        min(90, lat + target_cell_size_lat),
                    ),
                    (lon, min(90, lat + target_cell_size_lat)),
                    (lon, lat),
                ]
            )

            # Create record
            vgrid_record = graticule_dggs_to_geoseries(
                "vgrid", str(vgrid_id), resolution, cell_polygon
            )
            vgrid_records.append(vgrid_record)
            pbar.update(1)

    return gpd.GeoDataFrame(vgrid_records, geometry="geometry", crs="EPSG:4326")


def vgridgen(
    vgrid_instance: VGRID,
    resolution: int = 0,
    bbox: list[float] = None,
    output_format: str = "gpd",
):
    """
    Generate VGRID grid for pure Python usage.

    Note: Create a VGRID instance first using VGRID(cell_size, aperture).

    Args:
        vgrid_instance (VGRID): VGRID instance (create with VGRID(cell_size, aperture))
        resolution (int, optional): Resolution level to generate. If None, uses VGRID instance resolution.
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', 'gpd', 'shapefile', 'gpkg', 'parquet', or None for list of VGRID IDs).

    Returns:
        dict, list, or str: Output in the requested format or file path.
    """

    if bbox is None:
        # Calculate total cells for whole world
        total_columns = vgrid_instance.total_columns
        total_rows = vgrid_instance.total_rows
        total_cells = total_columns * total_rows

        if total_cells > MAX_CELLS:
            raise ValueError(
                f"VGRID instance will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
            )

    if output_format is None:
        # Return list of VGrid IDs
        return vgrid_gen_ids(vgrid_instance, resolution, bbox)
    else:
        # Return GeoDataFrame in specified format
        gdf = vgrid_gen(vgrid_instance, resolution, bbox)
        output_name = f"vgrid_grid_{resolution}"
        return convert_to_output_format(gdf, output_format, output_name)


def vgridgen_cli():
    """CLI interface for generating VGRID grid."""
    parser = argparse.ArgumentParser(description="Generate VGRID DGGS.")
    parser.add_argument(
        "-cell_size",
        "--cell_size",
        type=float,
        required=True,
        help="Cell size in degrees (must divide 360 and 180 evenly)",
    )
    parser.add_argument(
        "-aperture",
        "--aperture",
        type=int,
        default=4,
        help="Ratio to the next resolution (4 for quadtree-like, 9 for nonagon-like, default: 4)",
    )
    parser.add_argument(
        "-resolution",
        "--resolution",
        type=int,
        help="Resolution level to generate (optional, uses calculated resolution if not specified)",
    )
    parser.add_argument(
        "-bbox",
        "--bbox",
        nargs=4,
        type=float,
        metavar=("min_lon", "min_lat", "max_lon", "max_lat"),
        help="Bounding box coordinates (default: whole world)",
    )
    parser.add_argument(
        "-output_format",
        "--output_format",
        choices=OUTPUT_FORMATS,
        default="geojson",
        help="Output format (default: geojson)",
    )

    args = parser.parse_args()

    try:
        # Create VGRID instance
        vgrid_instance = VGRID(args.cell_size, args.aperture)

        result = vgridgen(
            vgrid_instance, args.resolution, args.bbox, args.output_format
        )
        if result is not None:
            print(f"VGRID grid generated successfully: {result}")
    except Exception as e:
        print(f"Error generating VGRID grid: {e}")
        return 1

    return 0


if __name__ == "__main__":
    vgridgen_cli()

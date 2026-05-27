"""
OLC Grid Generator Module

Generates OLC (Open Location Code) DGGS grids for specified resolutions with automatic cell generation and validation using human-readable location codes.

Key Functions:
- olc_grid(): Main grid generation function for whole world
- olc_grid_within_bbox(): Grid generation within bounding box
- olcgrid(): User-facing function with multiple output formats
- olcgrid_cli(): Command-line interface for grid generation
"""

import argparse
import geopandas as gpd
from vgrid.dggs import olc
from tqdm import tqdm
from shapely.geometry import shape, box, Polygon
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from shapely.ops import unary_union
from vgrid.utils.io import (
    is_full_world_bbox,
    validate_bbox,
    validate_olc_resolution,
    convert_to_output_format,
)


def olc_grid(resolution, verbose=True):
    resolution = validate_olc_resolution(resolution)
    """
    Generate a global grid of Open Location Codes (Plus Codes) at the specified precision
    as a GeoJSON-like feature collection.
    """
    # Define the boundaries of the world
    sw_lat, sw_lng = -90, -180
    ne_lat, ne_lng = 90, 180

    # Get the precision step size
    area = olc.decode(olc.encode(sw_lat, sw_lng, resolution))
    lat_step = area.latitudeHi - area.latitudeLo
    lng_step = area.longitudeHi - area.longitudeLo

    olc_records = []

    # Calculate the total number of steps for progress tracking
    total_lat_steps = int((ne_lat - sw_lat) / lat_step)
    total_lng_steps = int((ne_lng - sw_lng) / lng_step)
    total_steps = total_lat_steps * total_lng_steps

    with tqdm(
        total=total_steps,
        desc="Generating OLC DGGS",
        unit=" cells",
        disable=not verbose,
    ) as pbar:
        lat = sw_lat
        while lat < ne_lat:
            lng = sw_lng
            while lng < ne_lng:
                # Generate the Plus Code for the center of the cell
                center_lat = lat + lat_step / 2
                center_lon = lng + lng_step / 2
                olc_id = olc.encode(center_lat, center_lon, resolution)
                resolution = olc.decode(olc_id).codeLength
                cell_polygon = Polygon(
                    [
                        [lng, lat],  # SW
                        [lng, lat + lat_step],  # NW
                        [lng + lng_step, lat + lat_step],  # NE
                        [lng + lng_step, lat],  # SE
                        [lng, lat],  # Close the polygon
                    ]
                )
                olc_record = graticule_dggs_to_geoseries(
                    "olc", olc_id, resolution, cell_polygon
                )
                olc_records.append(olc_record)
                lng += lng_step
                pbar.update(1)  # Update progress bar
            lat += lat_step

    # Return the feature collection
    return gpd.GeoDataFrame(olc_records, geometry="geometry", crs="EPSG:4326")


def olc_grid_within_bbox(resolution, bbox):
    """
    Generate a grid of Open Location Codes (Plus Codes) within the specified bounding box.
    """
    resolution = validate_olc_resolution(resolution)
    bbox = validate_bbox(bbox)
    if is_full_world_bbox(bbox):
        return olc_grid(resolution)

    bbox_poly = box(*bbox)

    # Step 1: Generate base cells at the lowest resolution (e.g., resolution 2)
    base_resolution = 2
    base_gdf = olc_grid(base_resolution, verbose=False)

    # Step 2: Identify seed cells that intersect with the bounding box
    seed_cells = []
    for idx, base_cell in base_gdf.iterrows():
        base_cell_poly = base_cell["geometry"]
        if bbox_poly.intersects(base_cell_poly):
            seed_cells.append(base_cell)

    refined_records = []

    # Step 3: Iterate over seed cells and refine to the output resolution
    for seed_cell in seed_cells:
        seed_cell_poly = seed_cell["geometry"]

        if resolution == base_resolution:
            refined_records.append(seed_cell)
        else:
            # Refine the seed cell to the output resolution and add it to the output
            refined_records.extend(
                olc_refine_cell(
                    seed_cell_poly.bounds, base_resolution, resolution, bbox_poly
                )
            )

    gdf = gpd.GeoDataFrame(refined_records, geometry="geometry", crs="EPSG:4326")
    gdf = gdf[gdf["resolution"] == resolution]
    gdf = gdf.drop_duplicates(subset=["olc"])

    return gdf


def olc_refine_cell(bounds, current_resolution, target_resolution, bbox_poly):
    """
    Refine a cell defined by bounds to the target resolution, recursively refining intersecting cells.
    """
    min_lon, min_lat, max_lon, max_lat = bounds
    if current_resolution < 10:
        valid_resolution = current_resolution + 2
    else:
        valid_resolution = current_resolution + 1

    area = olc.decode(olc.encode(min_lat, min_lon, valid_resolution))
    lat_step = area.latitudeHi - area.latitudeLo
    lng_step = area.longitudeHi - area.longitudeLo

    olc_records = []
    lat = min_lat
    while lat < max_lat:
        lng = min_lon
        while lng < max_lon:
            # Define the bounds of the finer cell
            finer_cell_bounds = (lng, lat, lng + lng_step, lat + lat_step)
            finer_cell_poly = box(*finer_cell_bounds)

            if bbox_poly.intersects(finer_cell_poly):
                # Generate the Plus Code for the center of the finer cell
                center_lat = lat + lat_step / 2
                center_lon = lng + lng_step / 2
                olc_id = olc.encode(center_lat, center_lon, valid_resolution)
                resolution = olc.decode(olc_id).codeLength

                cell_polygon = Polygon(
                    [
                        [lng, lat],  # SW
                        [lng, lat + lat_step],  # NW
                        [lng + lng_step, lat + lat_step],  # NE
                        [lng + lng_step, lat],  # SE
                        [lng, lat],  # Close the polygon
                    ]
                )

                olc_record = graticule_dggs_to_geoseries(
                    "olc", olc_id, resolution, cell_polygon
                )
                olc_records.append(olc_record)

                # Recursively refine the cell if not at target resolution
                if valid_resolution < target_resolution:
                    olc_records.extend(
                        olc_refine_cell(
                            finer_cell_bounds,
                            valid_resolution,
                            target_resolution,
                            bbox_poly,
                        )
                    )

            lng += lng_step
            # pbar.update(1)
        lat += lat_step

    return olc_records


def olc_grid_ids(resolution):
    """
    Return a list of OLC (Plus Code) IDs for the whole world at the given resolution.
    """
    resolution = validate_olc_resolution(resolution)
    sw_lat, sw_lng = -90, -180
    ne_lat, ne_lng = 90, 180

    area = olc.decode(olc.encode(sw_lat, sw_lng, resolution))
    lat_step = area.latitudeHi - area.latitudeLo
    lng_step = area.longitudeHi - area.longitudeLo

    ids = []
    total_lat_steps = int((ne_lat - sw_lat) / lat_step)
    total_lng_steps = int((ne_lng - sw_lng) / lng_step)
    total_steps = total_lat_steps * total_lng_steps

    with tqdm(total=total_steps, desc="Generating OLC IDs", unit=" cells") as pbar:
        lat = sw_lat
        while lat < ne_lat:
            lng = sw_lng
            while lng < ne_lng:
                center_lat = lat + lat_step / 2
                center_lon = lng + lng_step / 2
                olc_id = olc.encode(center_lat, center_lon, resolution)
                ids.append(olc_id)
                lng += lng_step
                pbar.update(1)
            lat += lat_step

    return ids


def olc_grid_within_bbox_ids(resolution, bbox):
    """
    Return a list of OLC (Plus Code) IDs within a bounding box at the given resolution.
    """
    resolution = validate_olc_resolution(resolution)
    bbox = validate_bbox(bbox)
    if is_full_world_bbox(bbox):
        return olc_grid_ids(resolution)

    gdf = olc_grid_within_bbox(resolution, bbox)
    if gdf is None or gdf.empty:
        return []
    return list(gdf["olc"].drop_duplicates())


def olcgrid(resolution, bbox=None, output_format="gpd"):
    """
    Generate OLC grid for pure Python usage.

    Args:
        resolution (int): OLC resolution [2..15]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', 'geo', 'gpd', 'shapefile', 'gpkg', 'parquet', or None for list of OLC IDs).

    Returns:
        dict, list, or str: Output in the requested format or file path.
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        gdf = olc_grid(resolution)
    else:
        gdf = olc_grid_within_bbox(resolution, bbox)

    output_name = f"olc_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def olcgrid_cli():
    parser = argparse.ArgumentParser(description="Generate OLC DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [2..15]"
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
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    try:
        result = olcgrid(resolution, bbox, args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    olcgrid_cli()

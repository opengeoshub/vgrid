import argparse
import json
from shapely.geometry import Polygon, box
from tqdm import tqdm
from vgrid.utils.easedggs.constants import grid_spec, ease_crs, geo_crs, levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import (
    grid_ids_to_geos,
    geo_polygon_to_grid_ids,
)
from vgrid.generator.settings import max_cells, geodesic_dggs_to_feature

# Initialize the geodetic model

geo_bounds = grid_spec["geo"]
min_longitude = geo_bounds["min_x"]
min_lattitude = geo_bounds["min_y"]
max_longitude = geo_bounds["max_x"]
max_latitude = geo_bounds["max_y"]


def get_ease_cells(resolution):
    """
    Generate a list of cell IDs based on the resolution, row, and column.
    """
    n_row = levels_specs[resolution]["n_row"]
    n_col = levels_specs[resolution]["n_col"]

    # Generate list of cell IDs
    cell_ids = []

    # Loop through all rows and columns at the specified resolution
    for row in range(n_row):
        for col in range(n_col):
            # Generate base ID (e.g., L0.RRRCCC for res=0)
            base_id = f"L{resolution}.{row:03d}{col:03d}"

            # Add additional ".RC" for each higher resolution
            cell_id = base_id
            for i in range(1, resolution + 1):
                cell_id += f".{row:1d}{col:1d}"  # For res=1: L0.RRRCCC.RC, res=2: L0.RRRCCC.RC.RC, etc.

            # Append the generated cell ID to the list
            cell_ids.append(cell_id)

    return cell_ids


def get_ease_cells_bbox(resolution, bbox):
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt
    cells_bbox = geo_polygon_to_grid_ids(
        bounding_box_wkt,
        level=resolution,
        source_crs=geo_crs,
        target_crs=ease_crs,
        levels_specs=levels_specs,
        return_centroids=True,
        wkt_geom=True,
    )
    return cells_bbox


def generate_grid(resolution):
    ease_features = []

    level_spec = levels_specs[resolution]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]

    cells = get_ease_cells(resolution)

    for cell in tqdm(
        cells, total=len(cells), desc="Generating EASE DGGS", unit=" cells"
    ):
        geo = grid_ids_to_geos([cell])
        center_lon, center_lat = geo["result"]["data"][0]
        cell_min_lat = center_lat - (180 / (2 * n_row))
        cell_max_lat = center_lat + (180 / (2 * n_row))
        cell_min_lon = center_lon - (360 / (2 * n_col))
        cell_max_lon = center_lon + (360 / (2 * n_col))

        cell_polygon = Polygon(
            [
                [cell_min_lon, cell_min_lat],
                [cell_max_lon, cell_min_lat],
                [cell_max_lon, cell_max_lat],
                [cell_min_lon, cell_max_lat],
                [cell_min_lon, cell_min_lat],
            ]
        )
        if cell_polygon:
            num_edges = 4
            ease_feature = geodesic_dggs_to_feature(
                "ease", str(cell), resolution, cell_polygon, num_edges
            )
            ease_features.append(ease_feature)

    return {"type": "FeatureCollection", "features": ease_features}


def generate_grid_within_bbox(resolution, bbox):
    ease_features = []
    level_spec = levels_specs[resolution]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]

    # Get all grid cells within the bounding box
    cells = get_ease_cells_bbox(resolution, bbox)["result"]["data"]

    if cells:
        for cell in tqdm(cells, desc="Generating EASE DGGS", unit=" cells"):
            geo = grid_ids_to_geos([cell])
            if geo:
                center_lon, center_lat = geo["result"]["data"][0]
                cell_min_lat = center_lat - (180 / (2 * n_row))
                cell_max_lat = center_lat + (180 / (2 * n_row))
                cell_min_lon = center_lon - (360 / (2 * n_col))
                cell_max_lon = center_lon + (360 / (2 * n_col))

                cell_polygon = Polygon(
                    [
                        [cell_min_lon, cell_min_lat],
                        [cell_max_lon, cell_min_lat],
                        [cell_max_lon, cell_max_lat],
                        [cell_min_lon, cell_max_lat],
                        [cell_min_lon, cell_min_lat],
                    ]
                )
                num_edges = 4
                ease_feature = geodesic_dggs_to_feature(
                    "ease", str(cell), resolution, cell_polygon, num_edges
                )
                ease_features.append(ease_feature)

        return {"type": "FeatureCollection", "features": ease_features}


def main():
    parser = argparse.ArgumentParser(description="Generate EASE-DGGS DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="resolution [0..6]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)",
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = (
        args.bbox
        if args.bbox
        else [min_longitude, min_lattitude, max_longitude, max_latitude]
    )

    if resolution < 0 or resolution > 6:
        print("Please select a resolution in [0..6] range and try again")
        return

    if bbox == [min_longitude, min_lattitude, max_longitude, max_latitude]:
        level_spec = levels_specs[resolution]
        n_row = level_spec["n_row"]
        n_col = level_spec["n_col"]
        total_cells = n_row * n_col

        print(f"Resolution {resolution} will generate {total_cells} cells ")
        if total_cells > max_cells:
            print(f"which exceeds the limit of {max_cells}.")
            print("Please select a smaller resolution and try again.")
            return
        # Start generating and saving the grid in chunks
        geojson_features = generate_grid(resolution)

    else:
        geojson_features = generate_grid_within_bbox(resolution, bbox)

    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"ease_grid_{resolution}.geojson"
        with open(geojson_path, "w") as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()

import json
import argparse
from vgrid.utils import olc
from tqdm import tqdm
from shapely.geometry import box, Polygon
import locale
locale.setlocale(locale.LC_ALL, '')  # Use the system's default locale

def calculate_total_cells(code_length, bbox):
    """Calculate the total number of cells within the bounding box for a given resolution."""
    area = olc.decode(olc.encode(bbox[1], bbox[0], code_length))  # Use bbox min lat, min lon for the area
    lat_step = area.latitudeHi - area.latitudeLo
    lng_step = area.longitudeHi - area.longitudeLo

    sw_lng,sw_lat,ne_lng,ne_lat = bbox
    total_lat_steps = int((ne_lat - sw_lat) / lat_step)
    total_lng_steps = int((ne_lng - sw_lng) / lng_step)
    
    return total_lat_steps * total_lng_steps


def generate_grid(code_length):
    """
    Generate a global grid of Open Location Codes (Plus Codes) at the specified precision
    as a GeoJSON-like feature collection.
    """
    # Define the boundaries of the world
    sw_lat, sw_lng = -90, -180
    ne_lat, ne_lng = 90, 180

    # Get the precision step size
    area = olc.decode(olc.encode(sw_lat, sw_lng, code_length))
    lat_step = area.latitudeHi - area.latitudeLo
    lng_step = area.longitudeHi - area.longitudeLo

    features = []

    # Calculate the total number of steps for progress tracking
    total_lat_steps = int((ne_lat - sw_lat) / lat_step)
    total_lng_steps = int((ne_lng - sw_lng) / lng_step)
    total_steps = total_lat_steps * total_lng_steps

    # Iterate over the entire globe with tqdm for progress tracking
    with tqdm(total=total_steps, desc="Processing cells",unit=" cells") as pbar:
        lat = sw_lat
        while lat < ne_lat:
            lng = sw_lng
            while lng < ne_lng:
                # Generate the Plus Code for the center of the cell
                center_lat = lat + lat_step / 2
                center_lng = lng + lng_step / 2
                olc_code = olc.encode(center_lat, center_lng, code_length)
                resolution = olc.decode(olc_code).codeLength
                # Create the feature
                features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[
                            [lng, lat],  # SW
                            [lng, lat + lat_step],  # NW
                            [lng + lng_step, lat + lat_step],  # NE
                            [lng + lng_step, lat],  # SE
                            [lng, lat]  # Close the polygon
                        ]]
                    },
                    "properties": {
                        "olc": olc_code,
                        "resolution": resolution
                    }
                })

                lng += lng_step
                pbar.update(1)  # Update progress bar
            lat += lat_step

    # Return the feature collection
    return {
        "type": "FeatureCollection",
        "features": features
    }

def generate_grid_within_bbox(code_length, bbox):
    """
    Generate a grid of Open Location Codes (Plus Codes) within the specified bounding box.
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    bbox_poly = box(min_lon, min_lat, max_lon, max_lat)

    # Step 1: Generate base cells at the lowest resolution (e.g., resolution 2)
    base_resolution = 2
    base_cells = generate_grid(base_resolution)

    # Step 2: Identify seed cells that intersect with the bounding box
    seed_cells = []
    for base_cell in base_cells["features"]:
        base_cell_poly = Polygon(base_cell["geometry"]["coordinates"][0])
        if bbox_poly.intersects(base_cell_poly):
            seed_cells.append(base_cell)

    refined_features = []

    # Step 3: Iterate over seed cells and refine to the output resolution
    for seed_cell in seed_cells:
        seed_cell_poly = Polygon(seed_cell["geometry"]["coordinates"][0])

        if seed_cell_poly.contains(bbox_poly) and code_length == base_resolution:
            # Append the seed cell directly if fully contained and resolution matches
            refined_features.append(seed_cell)
        else:
            # Refine the seed cell to the output resolution and add it to the output
            refined_features.extend(
                refine_cell(seed_cell_poly.bounds, base_resolution, code_length, bbox_poly)
            )

    resolution_features = [
        feature for feature in refined_features if feature["properties"]["resolution"] == code_length
    ]

    final_features = []
    seen_olc_codes = set()  # Reset the set for final feature filtering

    for feature in resolution_features:
        olc_code = feature["properties"]["olc"]
        if olc_code not in seen_olc_codes:  # Check if OLC code is already in the set
            final_features.append(feature)
            seen_olc_codes.add(olc_code)


    return {
        "type": "FeatureCollection",
        "features": final_features
    }


def refine_cell(bounds, current_resolution, target_resolution, bbox_poly):
    """
    Refine a cell defined by bounds to the target resolution, recursively refining intersecting cells.
    """
    min_lon, min_lat, max_lon, max_lat = bounds
    if current_resolution < 10:
        valid_resolution = current_resolution + 2
    else: valid_resolution = current_resolution + 1
    
    area = olc.decode(olc.encode(min_lat, min_lon, valid_resolution))
    lat_step = area.latitudeHi - area.latitudeLo
    lng_step = area.longitudeHi - area.longitudeLo

    features = []
    total_lat_steps = int((max_lat - min_lat) / lat_step)
    total_lng_steps = int((max_lon - min_lon) / lng_step)
    total_steps = total_lat_steps * total_lng_steps

    with tqdm(total=total_steps, desc="Refining cells", unit=" cells") as pbar:
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
                    center_lng = lng + lng_step / 2
                    olc_code = olc.encode(center_lat, center_lng, valid_resolution)
                    resolution = olc.decode(olc_code).codeLength
                    
                    # Add the finer cell as a feature
                    features.append({
                        "type": "Feature",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[
                                [finer_cell_bounds[0], finer_cell_bounds[1]],  # SW
                                [finer_cell_bounds[0], finer_cell_bounds[3]],  # NW
                                [finer_cell_bounds[2], finer_cell_bounds[3]],  # NE
                                [finer_cell_bounds[2], finer_cell_bounds[1]],  # SE
                                [finer_cell_bounds[0], finer_cell_bounds[1]]  # Close the polygon
                            ]]
                        },
                        "properties": {
                            "olc": olc_code,
                            "resolution": resolution
                        }
                    })

                    # Recursively refine the cell if not at target resolution
                    if valid_resolution < target_resolution:
                        features.extend(
                            refine_cell(
                                finer_cell_bounds,
                                valid_resolution,
                                target_resolution,
                                bbox_poly
                            )
                        )

                lng += lng_step
                pbar.update(1)
            lat += lat_step

    return features

def validate_resolution(resolution):
    if resolution < 2 or resolution > 15:
        print(f"Please select a resolution in the range [2..15].")
        return False
    
    if resolution < 10 and resolution % 2 != 0:
        print(f"Please select an even resolution number if it is less than 10.")
        return False
    
    # If resolution is greater than or equal to 10, it can be either even or odd
    return True

def main():
    parser = argparse.ArgumentParser(description="Generate a grid of Open Location Codes (Plus Codes).")
    parser.add_argument(
        "-r", "--resolution",
        type=int,
        required=True,
        help="OLC code length/ resolution must be in [2..15] range, and ensure it is an even number when < 10."
    )
    
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    num_cells = calculate_total_cells(resolution, bbox)
    max_cells = 1_000_000
    if not validate_resolution(resolution):
        print(f"Please select a resolution in [2..15] range, and ensure it is an even number when < 10.")
        return

    if  bbox == [-180, -90, 180, 90]:
        print(f"The selected resolution will generate "
            f"{locale.format_string('%d', num_cells, grouping=True)} cells, ")

        if num_cells > max_cells:
            print(
                f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}."
            )
            print("Please select a smaller resolution and try again.")
            return
        geojson_features = generate_grid(resolution)
        geojson_path = f"olc_grid_{resolution}.geojson"
    else:
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        geojson_path = f"olc_grid_bbox_{resolution}.geojson"

    
    with open(geojson_path, "w") as f:
        json.dump(geojson_features, f, indent=2)
    
    print(f"OLC grid saved to {geojson_path}")

if __name__ == "__main__":
    main()

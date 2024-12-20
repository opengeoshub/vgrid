import json
import argparse
import vgrid.utils.olc as olc


def generate_grid(resolution, bounds):
    """
    Generate Plus Codes at a specific resolution and save as GeoJSON.

    Args:
        bounds (tuple): Bounding box as (min_lon, min_lat, max_lon, max_lat).
        resolution (int): Plus Code resolution (0 to 15).
    """
    min_lon, min_lat, max_lon, max_lat = bounds
    step = 20.0 / (10 ** (resolution // 2))  # Approximate step size for given resolution.

    geojson_features = []
    lat = min_lat
    while lat < max_lat:
        lon = min_lon
        while lon < max_lon:
            # Generate Plus Code
            olc_code = olc.encode(lat, lon, resolution)
            coord = olc.decode(olc_code)

            if coord:
                # Create bounding box coordinates for the polygon
                min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
                max_lat, max_lon = coord.latitudeHi, coord.longitudeHi

                # GeoJSON feature for the Plus Code
                feature = {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[
                            [min_lon, min_lat],
                            [min_lon, max_lat],
                            [max_lon, max_lat],
                            [max_lon, min_lat],
                            [min_lon, min_lat]
                        ]]
                    },
                    "properties": {
                        "plus_code": olc_code
                    }
                }
                geojson_features.append(feature)
            lon += step
        lat += step

    # Create GeoJSON structure
    return {
        "type": "FeatureCollection",
        "features": geojson_features
    }

 
def main():
    """
    Main function to parse command-line arguments and generate OLC grid.
    """
    parser = argparse.ArgumentParser(description="Generate Plus Codes (OLC) grid as GeoJSON.")
    parser.add_argument(
        '-r', '--resolution', type=int, required=True, 
        help="Resolution [0..15] of the grid."
    )
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)."
    )
    args = parser.parse_args()
    resolution = args.resolution
    # Set default bounding box to the whole world if not provided
    bounds = args.bbox if args.bbox else (-180, -90, 180, 90)

    # Generate the grid and save as GeoJSON
    geojson_features = generate_grid(resolution, bounds)
   # Save to GeoJSON file
    geojson_path = f"olc_grid_{resolution}.geojson"
    with open(geojson_path, "w") as f:
        json.dump(geojson_features, f, indent=2)
    
    print(f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()

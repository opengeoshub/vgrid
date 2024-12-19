# Reference: https://geohash.softeng.co/uekkn, https://github.com/vinsci/geohash, https://www.movable-type.co.uk/scripts/geohash.html?geohash=dp3
import  vgrid.utils.geohash as geohash
import argparse,locale,json
from shapely.geometry import Polygon, mapping
from tqdm import tqdm

locale.setlocale(locale.LC_ALL, '')
max_cells = 1000_000

def geohash_to_polygon(gh):
    """Convert geohash to a Shapely Polygon."""
    lat, lon = geohash.decode(gh)
    lat_err, lon_err = geohash.decode_exactly(gh)[2:]

    bbox = {
        'w': max(lon - lon_err, -180),
        'e': min(lon + lon_err, 180),
        's': max(lat - lat_err, -85.051129),
        'n': min(lat + lat_err, 85.051129)
    }

    return Polygon([
        (bbox['w'], bbox['s']),
        (bbox['w'], bbox['n']),
        (bbox['e'], bbox['n']),
        (bbox['e'], bbox['s']),
        (bbox['w'], bbox['s'])
    ])

def generate_grid(resolution):
    """Generate GeoJSON for the entire world at the given geohash resolution."""
    initial_geohashes = [
        "b", "c", "f", "g", "u", "v", "y", "z",
        "8", "9", "d", "e", "s", "t", "w", "x",
        "0", "1", "2", "3", "p", "q", "r", "k",
        "m", "n", "h", "j", "4", "5", "6", "7"
    ]

    def expand_geohash(gh, target_length, geohashes):
        if len(gh) == target_length:
            geohashes.add(gh)
            return
        for char in "0123456789bcdefghjkmnpqrstuvwxyz":
            expand_geohash(gh + char, target_length, geohashes)

    geohashes = set()
    for gh in initial_geohashes:
        expand_geohash(gh, resolution, geohashes)

    features = []
    for gh in tqdm(geohashes, desc="Generating grid", unit=" cells"):
        polygon = geohash_to_polygon(gh)
        features.append({
            "type": "Feature",
            "geometry": mapping(polygon),
            "properties": {
                "geohash": gh
            }
        })

    return {
        "type": "FeatureCollection",
        "features": features
    }


def generate_grid_within_bbox(resolution, bbox):
    """Generate GeoJSON for geohashes within a bounding box at the given resolution."""
    features = []

    # Step 1: Find the geohash covering the center of the bounding box
    bbox_center = ((bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2)
    center_geohash = geohash.encode(bbox_center[1], bbox_center[0], precision=resolution)

    # Step 2: Find the ancestor geohash that fully contains the bounding box
    def find_ancestor_geohash(center_geohash, bbox):
        for r in range(1, len(center_geohash) + 1):
            ancestor = center_geohash[:r]
            polygon = geohash_to_polygon(ancestor)
            if polygon.contains(Polygon.from_bounds(*bbox)):
                return ancestor
        return None  # Fallback if no ancestor is found

    ancestor_geohash = find_ancestor_geohash(center_geohash, bbox)

    if not ancestor_geohash:
        raise ValueError("No ancestor geohash fully contains the bounding box.")

    # Step 3: Expand geohashes recursively from the ancestor
    bbox_polygon = Polygon.from_bounds(*bbox)

    def expand_geohash(gh, target_length, geohashes):
        """Expand geohash only if it intersects the bounding box."""
        polygon = geohash_to_polygon(gh)
        if not polygon.intersects(bbox_polygon):
            return  # Skip this branch if it doesn't intersect the bounding box

        if len(gh) == target_length:
            geohashes.add(gh)  # Add to the set if it reaches the target resolution
            return

        for char in "0123456789bcdefghjkmnpqrstuvwxyz":
            expand_geohash(gh + char, target_length, geohashes)

    geohashes = set()
    expand_geohash(ancestor_geohash, resolution, geohashes)

    # Step 4: Generate features for geohashes that intersect the bounding box
    for gh in tqdm(geohashes, desc="Generating grid", unit=" cells"):
        polygon = geohash_to_polygon(gh)
        features.append({
            "type": "Feature",
            "geometry": mapping(polygon),
            "properties": {
                "geohash": gh
            }
        })

    return {
        "type": "FeatureCollection",
        "features": features
    }


def main():
    parser = argparse.ArgumentParser(description='Generate world polygons based on geohashes.')
    parser.add_argument(
        '-r', '--resolution', type=int, required=True,
        help='Resolution level for the geohashes (1-12)'
    )
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4,
        help='Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)'
    )
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox

    if not (1 <= resolution <= 12):
        print("Resolution must be between 1 and 12.")
        return

    # Validate resolution and calculate metrics
    if not bbox:
        num_cells = 32 ** resolution
        print(f"Resolution {resolution} will generate {locale.format_string('%d', num_cells, grouping=True)} cells.")
        if num_cells > max_cells:
            print(f"Exceeds limit of {locale.format_string('%d', max_cells, grouping=True)} cells.")
            print("Please select a smaller resolution and try again.")
            return

        geojson_features = generate_grid(resolution)

        # Define the GeoJSON file path
        geojson_path = f"geohash_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    
    else:
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"geohash_grid_{resolution}_bbox.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()



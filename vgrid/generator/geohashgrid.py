# Reference: https://geohash.softeng.co/uekkn, https://github.com/vinsci/geohash, https://www.movable-type.co.uk/scripts/geohash.html?geohash=dp3
import vgrid.utils.geohash as geohash
import argparse
import json
from shapely.geometry import Polygon, shape
from shapely.ops import unary_union
from tqdm import tqdm
from vgrid.generator.settings import max_cells, graticule_dggs_to_feature

initial_geohashes = [
    "b",
    "c",
    "f",
    "g",
    "u",
    "v",
    "y",
    "z",
    "8",
    "9",
    "d",
    "e",
    "s",
    "t",
    "w",
    "x",
    "0",
    "1",
    "2",
    "3",
    "p",
    "q",
    "r",
    "k",
    "m",
    "n",
    "h",
    "j",
    "4",
    "5",
    "6",
    "7",
]


def geohash_to_polygon(gh):
    """Convert geohash to a Shapely Polygon."""
    lat, lon = geohash.decode(gh)
    lat_err, lon_err = geohash.decode_exactly(gh)[2:]

    bbox = {
        "w": max(lon - lon_err, -180),
        "e": min(lon + lon_err, 180),
        "s": max(lat - lat_err, -85.051129),
        "n": min(lat + lat_err, 85.051129),
    }

    return Polygon(
        [
            (bbox["w"], bbox["s"]),
            (bbox["w"], bbox["n"]),
            (bbox["e"], bbox["n"]),
            (bbox["e"], bbox["s"]),
            (bbox["w"], bbox["s"]),
        ]
    )


def expand_geohash(gh, target_length, geohashes):
    if len(gh) == target_length:
        geohashes.add(gh)
        return
    for char in "0123456789bcdefghjkmnpqrstuvwxyz":
        expand_geohash(gh + char, target_length, geohashes)


def generate_grid(resolution):
    """Generate GeoJSON for the entire world at the given geohash resolution."""

    geohashes = set()
    for gh in initial_geohashes:
        expand_geohash(gh, resolution, geohashes)

    geohash_features = []
    for gh in tqdm(geohashes, desc="Generating Geohash DGGS", unit=" cells"):
        cell_polygon = geohash_to_polygon(gh)
        geohash_feature = graticule_dggs_to_feature(
            "geohash", gh, resolution, cell_polygon
        )
        geohash_features.append(geohash_feature)

    return {"type": "FeatureCollection", "features": geohash_features}


def expand_geohash_bbox(gh, target_length, geohashes, bbox_polygon):
    """Expand geohash only if it intersects the bounding box."""
    polygon = geohash_to_polygon(gh)
    if not polygon.intersects(bbox_polygon):
        return

    if len(gh) == target_length:
        geohashes.add(gh)  # Add to the set if it reaches the target resolution
        return

    for char in "0123456789bcdefghjkmnpqrstuvwxyz":
        expand_geohash_bbox(gh + char, target_length, geohashes, bbox_polygon)


def generate_grid_within_bbox(resolution, bbox):
    """Generate GeoJSON for geohashes within a bounding box at the given resolution."""
    geohash_features = []
    bbox_polygon = Polygon.from_bounds(*bbox)

    # Compute intersected geohashes using set comprehension
    intersected_geohashes = {
        gh
        for gh in initial_geohashes
        if geohash_to_polygon(gh).intersects(bbox_polygon)
    }

    # Expand geohash bounding box
    geohashes_bbox = set()
    for gh in intersected_geohashes:
        expand_geohash_bbox(gh, resolution, geohashes_bbox, bbox_polygon)

    # Generate GeoJSON features
    geohash_features.extend(
        graticule_dggs_to_feature("geohash", gh, resolution, geohash_to_polygon(gh))
        for gh in tqdm(geohashes_bbox, desc="Generating Geohash DGGS", unit=" cells")
    )

    return {"type": "FeatureCollection", "features": geohash_features}


def generate_grid_resample(resolution, geojson_features):
    """Generate GeoJSON for geohashes within a GeoJSON feature collection at the given resolution."""
    geohash_features = []

    # Union of all input geometries
    geometries = [
        shape(feature["geometry"]) for feature in geojson_features["features"]
    ]
    unified_geom = unary_union(geometries)

    # Compute intersected geohashes from the initial set
    intersected_geohashes = {
        gh
        for gh in initial_geohashes
        if geohash_to_polygon(gh).intersects(unified_geom)
    }

    # Expand geohash coverage within the unified geometry
    geohashes_geom = set()
    for gh in intersected_geohashes:
        expand_geohash_bbox(gh, resolution, geohashes_geom, unified_geom)

    # Generate GeoJSON features
    geohash_features.extend(
        graticule_dggs_to_feature("geohash", gh, resolution, geohash_to_polygon(gh))
        for gh in tqdm(geohashes_geom, desc="Generating Geohash DGGS", unit="cells")
    )

    return {"type": "FeatureCollection", "features": geohash_features}


def main():
    parser = argparse.ArgumentParser(description="Generate Geohash DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [1..10]"
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
    bbox = args.bbox

    if not (1 <= resolution <= 10):
        print("Resolution must be between 1 and 10.")
        return

    # Validate resolution and calculate metrics
    if not bbox:
        total_cells = 32**resolution
        print(f"Resolution {resolution} will generate {total_cells} cells ")
        if total_cells > max_cells:
            print(f"which exceeds the limit of {max_cells} cells.")
            print("Please select a smaller resolution and try again.")
            return

        geojson_features = generate_grid(resolution)

    else:
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)

    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"geohash_grid_{resolution}.geojson"
        with open(geojson_path, "w") as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()

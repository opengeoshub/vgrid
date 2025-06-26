from vgrid.utils import geohash
from shapely.geometry import Point, LineString, Polygon
import argparse
import json
from tqdm import tqdm
import os
import requests
from urllib.parse import urlparse
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.generator.geohashgrid import (
    initial_geohashes,
    geohash_to_polygon,
    expand_geohash_bbox,
)
from vgrid.conversion.dggscompact import geohashcompact


# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties):
    geohash_features = []
    longitude = point.x
    latitude = point.y
    geohash_id = geohash.encode(latitude, longitude, resolution)
    bbox = geohash.bbox(geohash_id)
    if bbox:
        min_lat, min_lon = bbox["s"], bbox["w"]  # Southwest corner
        max_lat, max_lon = bbox["n"], bbox["e"]  # Northeast corner
        resolution = len(geohash_id)

        # Define the polygon based on the bounding box
        cell_polygon = Polygon(
            [
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat],  # Closing the polygon (same as the first point)
            ]
        )
        geohash_feature = graticule_dggs_to_feature(
            "geohash", geohash_id, resolution, cell_polygon
        )
        geohash_feature["properties"].update(feature_properties)
        geohash_features.append(geohash_feature)

    return {"type": "FeatureCollection", "features": geohash_features}


def poly_to_grid(resolution, geometry, feature_properties, compact):
    geohash_features = []
    if geometry.geom_type == "LineString" or geometry.geom_type == "Polygon":
        polys = [geometry]
    elif (
        geometry.geom_type == "MultiLineString" or geometry.geom_type == "MultiPolygon"
    ):
        polys = list(geometry)

    for poly in polys:
        intersected_geohashes = {
            gh for gh in initial_geohashes if geohash_to_polygon(gh).intersects(poly)
        }
        # Expand geohash bounding box
        geohashes_bbox = set()
        for gh in intersected_geohashes:
            expand_geohash_bbox(gh, resolution, geohashes_bbox, poly)

        # Process geohashes
        for gh in geohashes_bbox:
            cell_polygon = geohash_to_polygon(gh)
            geohash_feature = graticule_dggs_to_feature(
                "geohash", gh, resolution, cell_polygon
            )
            geohash_feature["properties"].update(feature_properties)

            geohash_features.append(geohash_feature)

    geohash_geosjon = {"type": "FeatureCollection", "features": geohash_features}

    if compact:
        return geohashcompact(geohash_geosjon)

    return geohash_geosjon


def geojson2geohash(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to Geohash DGGS format.

    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): Resolution level [1..10]
        compact (bool): Whether to enable Geohash compact mode

    Returns:
        dict: GeoJSON data with Geohash features
    """
    if resolution < 1 or resolution > 10:
        raise ValueError("Resolution must be in range [1..10]")

    geojson_features = []

    for feature in tqdm(geojson_data["features"], desc="Processing GeoJSON features"):
        feature_properties = feature["properties"]
        if feature["geometry"]["type"] in ["Point", "MultiPoint"]:
            coordinates = feature["geometry"]["coordinates"]
            if feature["geometry"]["type"] == "Point":
                point = Point(coordinates)
                point_features = point_to_grid(resolution, point, feature_properties)
                geojson_features.extend(point_features["features"])

            elif feature["geometry"]["type"] == "MultiPoint":
                for point_coords in coordinates:
                    point = Point(point_coords)
                    point_features = point_to_grid(
                        resolution, point, feature_properties
                    )
                    geojson_features.extend(point_features["features"])

        elif feature["geometry"]["type"] in ["LineString", "MultiLineString"]:
            coordinates = feature["geometry"]["coordinates"]
            if feature["geometry"]["type"] == "LineString":
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(
                    resolution, polyline, feature_properties, compact
                )
                geojson_features.extend(polyline_features["features"])

            elif feature["geometry"]["type"] == "MultiLineString":
                for line_coords in coordinates:
                    polyline = LineString(line_coords)
                    polyline_features = poly_to_grid(
                        resolution, polyline, feature_properties, compact
                    )
                    geojson_features.extend(polyline_features["features"])

        elif feature["geometry"]["type"] in ["Polygon", "MultiPolygon"]:
            coordinates = feature["geometry"]["coordinates"]

            if feature["geometry"]["type"] == "Polygon":
                exterior_ring = coordinates[0]
                interior_rings = coordinates[1:]
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(
                    resolution, polygon, feature_properties, compact
                )
                geojson_features.extend(polygon_features["features"])

            elif feature["geometry"]["type"] == "MultiPolygon":
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(
                        resolution, polygon, feature_properties, compact
                    )
                    geojson_features.extend(polygon_features["features"])

    return {"type": "FeatureCollection", "features": geojson_features}


def is_url(path):
    """Check if the given path is a URL."""
    try:
        result = urlparse(path)
        return all([result.scheme, result.netloc])
    except Exception:
        return False


def read_geojson_file(geojson_path):
    """Read GeoJSON from either a local file or URL."""
    if is_url(geojson_path):
        try:
            response = requests.get(geojson_path)
            response.raise_for_status()
            return json.loads(response.text)
        except requests.RequestException as e:
            print(
                f"Error: Failed to download GeoJSON from URL {geojson_path}: {str(e)}"
            )
            return None
    else:
        if not os.path.exists(geojson_path):
            print(f"Error: The file {geojson_path} does not exist.")
            return None
        try:
            with open(geojson_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"Error reading GeoJSON file: {e}")
            return None


def geojson2geohash_cli():
    """
    Command-line interface for converting GeoJSON to Geohash DGGS format.
    Supports both local files and remote URLs.
    """
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Geohash DGGS")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [1..10]"
    )
    parser.add_argument(
        "-geojson",
        "--geojson",
        type=str,
        required=True,
        help="GeoJSON file path or URL (Point, Polyline or Polygon)",
    )
    parser.add_argument(
        "-compact", action="store_true", help="Enable Geohash compact mode"
    )

    args = parser.parse_args()

    # Read GeoJSON data from file or URL
    geojson_data = read_geojson_file(args.geojson)
    if geojson_data is None:
        return

    try:
        result = geojson2geohash(geojson_data, args.resolution, args.compact)

        # Save the results to GeoJSON
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2geohash_{args.resolution}.geojson"
        if args.compact:
            geojson_path = f"{geojson_name}2geohash_{args.resolution}_compacted.geojson"

        with open(geojson_path, "w") as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")

    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

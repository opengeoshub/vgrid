from shapely.geometry import Point, LineString, Polygon
import argparse
import json
from tqdm import tqdm
import os
import requests
from urllib.parse import urlparse
from vgrid.conversion.latlon2dggs import latlon2mgrs
from vgrid.conversion.dggs2geojson import mgrs2geojson


def point_to_grid(resolution, point, feature_properties):
    mgrs_features = []
    latitude, longitude = point.y, point.x
    mgrs_id = latlon2mgrs(latitude, longitude, resolution)
    mgrs_feature = mgrs2geojson(mgrs_id)

    if mgrs_feature:
        mgrs_feature["properties"].update(feature_properties)
        mgrs_features.append(mgrs_feature)

    return {
        "type": "FeatureCollection",
        "features": mgrs_features,
    }


# Function to generate grid for Polyline
def poly_to_grid(resolution, geometry, feature_properties):
    mgrs_features = []
    return {"type": "FeatureCollection", "features": mgrs_features}


def geojson2mgrs(geojson_data, resolution):
    """
    Convert GeoJSON data to MGRS DGGS format.

    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): MGRS resolution [0..5]

    Returns:
        dict: GeoJSON FeatureCollection with MGRS grid cells
    """
    if resolution < 0 or resolution > 5:
        raise ValueError("Resolution must be in range [0..5]")

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
                    point = Point(point_coords)  # Create Point for each coordinate set
                    point_features = point_to_grid(
                        resolution, point, feature_properties
                    )
                    geojson_features.extend(point_features["features"])

        elif feature["geometry"]["type"] in ["LineString", "MultiLineString"]:
            coordinates = feature["geometry"]["coordinates"]
            if feature["geometry"]["type"] == "LineString":
                # Directly process LineString geometry
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(
                    resolution, polyline, feature_properties
                )
                geojson_features.extend(polyline_features["features"])

            elif feature["geometry"]["type"] == "MultiLineString":
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = poly_to_grid(
                        resolution, polyline, feature_properties
                    )
                    geojson_features.extend(polyline_features["features"])

        elif feature["geometry"]["type"] in ["Polygon", "MultiPolygon"]:
            coordinates = feature["geometry"]["coordinates"]

            if feature["geometry"]["type"] == "Polygon":
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[
                    0
                ]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[
                    1:
                ]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(resolution, polygon, feature_properties)
                geojson_features.extend(polygon_features["features"])

            elif feature["geometry"]["type"] == "MultiPolygon":
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[
                        0
                    ]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[
                        1:
                    ]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(
                        resolution, polygon, feature_properties
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


def geojson2mgrs_cli():
    """
    Command-line interface for converting GeoJSON to MGRS DGGS format.
    Supports both local files and remote URLs.
    """
    parser = argparse.ArgumentParser(description="Convert GeoJSON to MGRS DGGS")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..5]"
    )
    parser.add_argument(
        "-geojson",
        "--geojson",
        type=str,
        required=True,
        help="GeoJSON file path or URL (Point, Polyline or Polygon)",
    )
    args = parser.parse_args()

    # Read GeoJSON data from file or URL
    geojson_data = read_geojson_file(args.geojson)
    if geojson_data is None:
        return

    try:
        result = geojson2mgrs(geojson_data, args.resolution)

        # Save the results to GeoJSON
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2mgrs_{args.resolution}.geojson"
        with open(geojson_path, "w") as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")

    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

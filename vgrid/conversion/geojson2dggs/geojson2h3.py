import os
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon, box, MultiPoint
import h3
import requests
from urllib.parse import urlparse
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.geometry import calculate_point_distances
from pyproj import Geod

geod = Geod(ellps="WGS84")


def get_nearest_h3_resolution(points):
    """
    Find the first H3 resolution where avg_edge_length < shortest_distance.
    Uses h3.average_edge_length(resolution) to compare with the calculated shortest distance.

    Args:
        points: Shapely Point geometry or MultiPoint geometry

    Returns:
        int: First H3 resolution where avg_edge_length < shortest_distance (0-15)
    """
    # Calculate the shortest distance between points
    shortest_distance = calculate_point_distances(points)

    # If only one point or no distance, return a default resolution
    if shortest_distance == 0:
        return 0  # Default resolution for single points

    # Check resolutions from 0 to 15 (lowest to highest)
    for resolution in range(16):
        # Get average edge length for this H3 resolution in meters
        avg_edge_length = h3.average_hexagon_edge_length(res=resolution, unit="m")

        # Return the first resolution where avg_edge_length < shortest_distance
        if avg_edge_length <= shortest_distance:
            return resolution

    # If no resolution found where avg_edge_length < shortest_distance, return the highest resolution
    return 15


# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties, topology=False):
    h3_features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    h3_id = h3.latlng_to_cell(latitude, longitude, resolution)

    cell_boundary = h3.cell_to_boundary(h3_id)
    # Wrap and filter the boundary
    filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
    # Reverse lat/lon to lon/lat for GeoJSON compatibility
    reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
    cell_polygon = Polygon(reversed_boundary)
    if cell_polygon:
        num_edges = 6
        if h3.is_pentagon(h3_id):
            num_edges = 5
        h3_feature = geodesic_dggs_to_feature(
            "h3", h3_id, resolution, cell_polygon, num_edges
        )
        h3_feature["properties"].update(feature_properties)
        h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def geodesic_buffer(polygon, distance):
    buffered_coords = []
    for lon, lat in polygon.exterior.coords:
        # Generate points around the current vertex to approximate a circle
        circle_coords = [
            geod.fwd(lon, lat, azimuth, distance)[
                :2
            ]  # Forward calculation: returns (lon, lat, back_azimuth)
            for azimuth in range(0, 360, 10)  # Generate points every 10 degrees
        ]
        buffered_coords.append(circle_coords)

    # Flatten the list of buffered points and form a Polygon
    all_coords = [coord for circle in buffered_coords for coord in circle]
    return Polygon(all_coords).convex_hull


def poly_to_grid(
    resolution, geometry, feature_properties, compact=False, topology=False
):
    h3_features = []

    if geometry.geom_type == "LineString" or geometry.geom_type == "Polygon":
        polys = [geometry]
    elif (
        geometry.geom_type == "MultiLineString" or geometry.geom_type == "MultiPolygon"
    ):
        polys = list(geometry)

    for poly in polys:
        bbox = box(*poly.bounds)
        distance = h3.average_hexagon_edge_length(resolution, unit="m") * 2
        bbox_buffer = geodesic_buffer(bbox, distance)
        bbox_buffer_cells = h3.geo_to_cells(bbox_buffer, resolution)
        if compact:
            bbox_buffer_cells = h3.compact_cells(bbox_buffer_cells)

        for bbox_buffer_cell in bbox_buffer_cells:
            cell_boundary = h3.cell_to_boundary(bbox_buffer_cell)
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            if cell_polygon.intersects(poly):
                cell_resolution = h3.get_resolution(bbox_buffer_cell)
                num_edges = 6
                if h3.is_pentagon(bbox_buffer_cell):
                    num_edges = 5
                h3_feature = geodesic_dggs_to_feature(
                    "h3", bbox_buffer_cell, cell_resolution, cell_polygon, num_edges
                )
                h3_feature["properties"].update(feature_properties)
                h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def geojson2h3(geojson_data, resolution, compact=False, topology=False):
    """
    Convert GeoJSON data to H3 grid cells.

    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): H3 resolution [0..15]
        compact (bool): Enable H3 compact mode - for polygon only
        topology (bool): Enable H3 topology preserving mode

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    if resolution < 0 or resolution > 15:
        raise ValueError("Resolution must be in range [0..15]")

    geojson_features = []

    # If topology mode is enabled, convert GeoJSON to Shapely geometries first
    if topology:
        # Extract all geometries from GeoJSON
        geometries = []
        for feature in geojson_data["features"]:
            geom_type = feature["geometry"]["type"]
            coords = feature["geometry"]["coordinates"]

            if geom_type == "Point":
                geometries.append(Point(coords))
            elif geom_type == "MultiPoint":
                for point_coords in coords:
                    geometries.append(Point(point_coords))
            elif geom_type == "LineString":
                geometries.append(LineString(coords))
            elif geom_type == "MultiLineString":
                for line_coords in coords:
                    geometries.append(LineString(line_coords))
            elif geom_type == "Polygon":
                exterior_ring = coords[0]
                interior_rings = coords[1:]
                geometries.append(Polygon(exterior_ring, interior_rings))
            elif geom_type == "MultiPolygon":
                for sub_polygon_coords in coords:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    geometries.append(Polygon(exterior_ring, interior_rings))

        # Create MultiPoint from all geometries for resolution calculation
        all_points = []
        for geom in geometries:
            if geom.geom_type == "Point":
                all_points.append(geom)
            elif geom.geom_type == "MultiPoint":
                all_points.extend(list(geom.geoms))
            else:
                # For other geometry types, add representative points
                all_points.append(geom.representative_point())

        if all_points:
            points = MultiPoint(all_points)
            resolution = get_nearest_h3_resolution(points)

    for feature in tqdm(geojson_data["features"], desc="Processing GeoJSON features"):
        feature_properties = feature["properties"]
        if feature["geometry"]["type"] in ["Point", "MultiPoint"]:
            coordinates = feature["geometry"]["coordinates"]
            if feature["geometry"]["type"] == "Point":
                point = Point(coordinates)
                point_features = point_to_grid(
                    resolution, point, feature_properties, topology
                )
                geojson_features.extend(point_features["features"])

            elif feature["geometry"]["type"] == "MultiPoint":
                for point_coords in coordinates:
                    point = Point(point_coords)
                    point_features = point_to_grid(
                        resolution, point, feature_properties, topology
                    )
                    geojson_features.extend(point_features["features"])

        elif feature["geometry"]["type"] in ["LineString", "MultiLineString"]:
            coordinates = feature["geometry"]["coordinates"]
            if feature["geometry"]["type"] == "LineString":
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(
                    resolution, polyline, feature_properties, topology
                )
                geojson_features.extend(polyline_features["features"])

            elif feature["geometry"]["type"] == "MultiLineString":
                for line_coords in coordinates:
                    polyline = LineString(line_coords)
                    polyline_features = poly_to_grid(
                        resolution, polyline, feature_properties
                    )
                    geojson_features.extend(polyline_features["features"])

        elif feature["geometry"]["type"] in ["Polygon", "MultiPolygon"]:
            coordinates = feature["geometry"]["coordinates"]

            if feature["geometry"]["type"] == "Polygon":
                exterior_ring = coordinates[0]
                interior_rings = coordinates[1:]
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(
                    resolution, polygon, feature_properties, compact, topology
                )
                geojson_features.extend(polygon_features["features"])

            elif feature["geometry"]["type"] == "MultiPolygon":
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(
                        resolution, polygon, feature_properties, compact, topology
                    )
                    geojson_features.extend(polygon_features["features"])

    return {
        "type": "FeatureCollection",
        "features": geojson_features,
    }


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


def geojson2h3_cli():
    """
    Command-line interface for converting GeoJSON to H3 grid cells.
    Supports both local files and remote URLs.
    """
    parser = argparse.ArgumentParser(description="Convert GeoJSON to H3 DGGS")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..15]"
    )
    parser.add_argument(
        "-geojson",
        "--geojson",
        type=str,
        required=True,
        help="GeoJSON file path or URL (Point, Polyline or Polygon)",
    )
    parser.add_argument(
        "-compact",
        action="store_true",
        help="Enable H3 compact mode - for polygon only",
    )
    parser.add_argument(
        "-topology", action="store_true", help="Enable H3 topology preserving mode"
    )

    args = parser.parse_args()

    # Read GeoJSON data from file or URL
    geojson_data = read_geojson_file(args.geojson)
    if geojson_data is None:
        return

    try:
        result = geojson2h3(geojson_data, args.resolution, args.compact, args.topology)

        # Save the results to GeoJSON
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2h3_{args.resolution}.geojson"
        suffix = ""
        if args.compact:
            suffix = "_compacted"
        if args.topology:
            suffix += "_topology"
        geojson_path = f"{geojson_name}2h3_{args.resolution}{suffix}.geojson"

        with open(geojson_path, "w") as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")

    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

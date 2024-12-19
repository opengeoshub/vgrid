import json
import numpy as np
from scipy.spatial import SphericalVoronoi
from scipy.spatial.transform import Rotation as R
from shapely.geometry import shape, Polygon, MultiPolygon
from vgrid.utils.antimeridian import fix_polygon  # Assumes fix_polygon is implemented

def fix_antimeridian_cells(hex_boundary, threshold=-128):
    """
    Adjust hex boundary longitudes to handle antimeridian crossings.
    """
    if any(lon < threshold for _, lon in hex_boundary):
        # Adjust all longitudes accordingly
        return [(lat, lon - 360 if lon > 0 else lon) for lat, lon in hex_boundary]
    return hex_boundary


def geojson_to_cartesian(points, radius=1):
    """
    Convert geographic coordinates to Cartesian coordinates on a sphere.
    Latitude and Longitude are expected in degrees.
    """
    lat = np.radians(points[:, 1])  # Latitude
    lon = np.radians(points[:, 0])  # Longitude
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.column_stack((x, y, z))


def cartesian_to_geojson(vertices, radius=1):
    """
    Convert Cartesian coordinates on a sphere to geographic coordinates.
    """
    vertices = np.array(vertices)
    lat = np.arcsin(vertices[:, 2] / radius)
    lon = np.arctan2(vertices[:, 1], vertices[:, 0])
    return np.degrees(lon), np.degrees(lat)


def extract_polygon_centroids(geojson_input):
    """
    Extract centroids of polygons from a GeoJSON file.
    """
    with open(geojson_input) as f:
        geojson_data = json.load(f)

    centroids = []
    for feature in geojson_data["features"]:
        geom = shape(feature["geometry"])
        if geom.geom_type in {"Polygon", "MultiPolygon"}:
            centroid = geom.centroid
            centroids.append([centroid.x, centroid.y])

    return np.array(centroids)  # Return as a NumPy array for easier processing


def create_spherical_voronoi(geojson_input, radius=1):
    """
    Create a spherical Voronoi diagram from GeoJSON polygon centroids,
    fixing antimeridian issues.
    """
    # Step 1: Extract centroids from GeoJSON polygons
    points = extract_polygon_centroids(geojson_input)

    # Step 2: Convert points to Cartesian coordinates
    cartesian_points = geojson_to_cartesian(points, radius=radius)

    # Step 3: Generate Voronoi diagram
    sv = SphericalVoronoi(cartesian_points, radius=radius)
    sv.sort_vertices_of_regions()

    # Step 4: Convert Voronoi vertices back to geographic coordinates and fix antimeridian issues
    voronoi_geojson = {
        "type": "FeatureCollection",
        "features": []
    }

    for region in sv.regions:
        vertices = sv.vertices[region]
        lon, lat = cartesian_to_geojson(vertices, radius)
        hex_boundary = list(zip(lat, lon))  # Create a list of (lat, lon) tuples

        # Fix antimeridian crossings
        fixed_boundary = fix_antimeridian_cells(hex_boundary)

        # Convert fixed boundary back to GeoJSON-compatible coordinates
        fixed_coordinates = [[lon, lat] for lat, lon in fixed_boundary]
        fixed_coordinates.append(fixed_coordinates[0])  # Close the polygon

        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [fixed_coordinates]
            },
            "properties": {}
        }
        voronoi_geojson["features"].append(feature)

    return voronoi_geojson


# Example usage
if __name__ == "__main__":
    geojson_input = "rhealpix_grid_2.geojson"  # Replace with your GeoJSON file containing polygons
    radius = 6371  # Radius of the Earth in kilometers
    voronoi_result = create_spherical_voronoi(geojson_input, radius)

    # Save result to a new GeoJSON file
    with open("voronoi.geojson", "w") as f:
        json.dump(voronoi_result, f)

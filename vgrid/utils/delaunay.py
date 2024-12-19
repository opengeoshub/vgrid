import json
import numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon, shape, mapping
from shapely.ops import unary_union
from typing import Union
from vgrid.utils.antimeridian import fix_polygon

def spherical_to_cartesian(lon, lat, radius=1):
    """
    Convert spherical coordinates to Cartesian coordinates.
    Longitude and latitude are in degrees.
    """
    lon = np.radians(lon)
    lat = np.radians(lat)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.stack((x, y, z), axis=-1)


def cartesian_to_spherical(x, y, z):
    """
    Convert Cartesian coordinates to spherical coordinates.
    Returns longitude and latitude in degrees.
    """
    lon = np.degrees(np.arctan2(y, x))
    lat = np.degrees(np.arcsin(z / np.sqrt(x**2 + y**2 + z**2)))
    return lon, lat


def extract_polygon_centroids(geojson_file):
    """
    Extract centroids of polygons from a GeoJSON file.
    """
    with open(geojson_file) as f:
        geojson_data = json.load(f)

    centroids = []
    for feature in geojson_data["features"]:
        geom = shape(feature["geometry"])
        if geom.geom_type in {"Polygon", "MultiPolygon"}:
            centroid = geom.centroid
            centroids.append([centroid.x, centroid.y])

    return centroids


def create_spherical_delaunay_polygons_from_geojson(geojson_file, radius=1):
    """
    Create a spherical Delaunay triangulation using centroids of polygons from a GeoJSON file,
    and return polygons with antimeridian adjustments.

    Parameters:
        geojson_file: Path to the GeoJSON file containing polygons.
        radius: Radius of the sphere (default is 1).

    Returns:
        List of Delaunay triangles as Polygons, adjusted for antimeridian crossings.
    """
    # Step 1: Extract centroids from GeoJSON
    centroids = extract_polygon_centroids(geojson_file)

    # Step 2: Convert centroids to Cartesian coordinates
    cartesian_points = spherical_to_cartesian(
        lon=[point[0] for point in centroids],
        lat=[point[1] for point in centroids],
        radius=radius
    )

    # Step 3: Perform Delaunay triangulation in Cartesian coordinates
    delaunay = Delaunay(cartesian_points)

    # Step 4: Extract triangles from the Delaunay triangulation
    polygons = []
    for simplex in delaunay.simplices:
        triangle_cartesian = cartesian_points[simplex]
        lons, lats = zip(*[cartesian_to_spherical(*point) for point in triangle_cartesian])

        # Create a Polygon
        triangle_polygon = Polygon(zip(lons, lats))

        # Fix antimeridian issues in the Polygon
        fixed_polygon = fix_polygon(triangle_polygon)
        polygons.append(triangle_polygon)

    return polygons


# Example usage
if __name__ == "__main__":
    geojson_file = "rhealpix_grid_2.geojson"  # Path to your GeoJSON file
    radius = 6371  # Earth's radius in kilometers
    delaunay_polygons = create_spherical_delaunay_polygons_from_geojson(geojson_file, radius)

    # Save the Polygons to a new GeoJSON file
    output_geojson = {
        "type": "FeatureCollection",
        "features": []
    }
    for polygon in delaunay_polygons:
        feature = {
            "type": "Feature",
            "geometry": mapping(polygon),
            "properties": {}
        }
        output_geojson["features"].append(feature)

    with open("delaunay.geojson", "w") as f:
        json.dump(output_geojson, f)

import json
import os
from shapely.geometry import Polygon, MultiPolygon, mapping
from shapely.validation import make_valid
import math
from vgrid.utils.antimeridian import fix_polygon

# # Function to fix invalid polygons
# def fix_polygon(polygon):
#     if not polygon.is_valid:
#         return make_valid(polygon)
#     return polygon

# Function to construct geometry
def construct_geometry(coords, is_multipolygon=False):
    if is_multipolygon:
        polygon = MultiPolygon([Polygon(poly) for poly in coords])
    else:
        polygon = Polygon(coords)

    fixed_polygon = fix_polygon(polygon)    
    return fixed_polygon
    # return polygon

# Convert 3D Cartesian coordinates to WGS84 (latitude, longitude)
def cartesian_to_wgs84(x, y, z):
    lon = math.degrees(math.atan2(y, x))
    hyp = math.hypot(x, y)
    lat = math.degrees(math.atan2(z, hyp))
    return (lon, lat)

# Function to generate global icosahedron
def main():
    out_file_dir = os.getcwd()
    out_file_name = "icosahedron.geojson"
    out_file = os.path.join(out_file_dir, out_file_name)

    # Golden ratio for icosahedron coordinates
    phi = (1 + 5 ** 0.5) / 2

    # Define vertices of icosahedron in 3D Cartesian coordinates
    vertices = [
        (0, 1, phi), (0, -1, phi), (0, 1, -phi), (0, -1, -phi),
        (1, phi, 0), (-1, phi, 0), (1, -phi, 0), (-1, -phi, 0),
        (phi, 0, 1), (-phi, 0, 1), (phi, 0, -1), (-phi, 0, -1)
    ]

    # Normalize vertices to unit sphere and convert to WGS84
    vertices = [cartesian_to_wgs84(x, y, z) for x, y, z in vertices]

    # Define faces by vertex indices
    faces_indices = [
        (0, 1, 8), (0, 4, 5), (0, 5, 9), (0, 8, 4), (0, 9, 1),
        (1, 6, 8), (1, 7, 6), (1, 9, 7), (2, 3, 10), (2, 4, 8),
        (2, 5, 4), (2, 10, 5), (3, 6, 7), (3, 7, 11), (3, 11, 10),
        (4, 8, 6), (4, 6, 5), (5, 6, 9), (7, 9, 11), (9, 10, 11)
    ]

    # Create GeoJSON features
    geojson_features = []
    for idx, (i, j, k) in enumerate(faces_indices, start=1):
        coords = [vertices[i], vertices[j], vertices[k], vertices[i]]
        geometry = construct_geometry(coords)
        geojson_features.append({
            "type": "Feature",
            "geometry": mapping(geometry),
            "properties": {"zoneID": str(idx)}
        })

    # Final GeoJSON output
    geojson_output = {
        "type": "FeatureCollection",
        "features": geojson_features
    }

    # Save as GeoJSON
    with open(out_file, 'w') as f:
        json.dump(geojson_output, f, indent=2)

    print(f"Icosahedron saved as {out_file}")

if __name__ == '__main__':
    main()

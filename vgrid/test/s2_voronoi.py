
import vgrid.utils.s2 as s2sphere
import geojson
from shapely.geometry import Polygon, mapping

# Define generating points (latitude, longitude pairs)
points = [
    s2sphere.LatLng.from_degrees(40.7128, -74.0060),  # New York
    s2sphere.LatLng.from_degrees(34.0522, -118.2437), # Los Angeles
    s2sphere.LatLng.from_degrees(51.5074, -0.1278)    # London
]

# Create an S2RegionCoverer to compute Voronoi cells
region_coverer = s2sphere.RegionCoverer()
region_coverer.min_level = 6  # Set the resolution (lower = coarser, higher = finer)
region_coverer.max_level = 6  # Keep the same level for consistent tiles
region_coverer.max_cells = 1000000

# Generate Voronoi cells
features = []
for i, point in enumerate(points):
    # Create an S2Cap around the point to simulate the Voronoi cell
    cap = s2sphere.Cap.from_axis_angle(point.to_point(), s2sphere.Angle.from_degrees(20))  # Adjustable angle
    covering = region_coverer.get_covering(cap)
    
    # Convert S2CellUnion to polygons
    for cell_id in covering:
        cell = s2sphere.Cell(cell_id)
        vertices = []
        for j in range(4):  # Each S2 cell has 4 vertices
            vertex = cell.get_vertex(j)
            lat_lng = s2sphere.LatLng.from_point(vertex)
            vertices.append((lat_lng.lng().degrees, lat_lng.lat().degrees))  # (lon, lat)
        vertices.append(vertices[0])  # Close the polygon
        
        # Create a GeoJSON feature for the cell
        polygon = Polygon(vertices)
        features.append(geojson.Feature(geometry=mapping(polygon), properties={"point_id": i}))

# Create a GeoJSON FeatureCollection
feature_collection = geojson.FeatureCollection(features)

# Save to GeoJSON file
with open("voronoi_sphere_s2.geojson", "w") as f:
    geojson.dump(feature_collection, f, indent=2)

print("GeoJSON saved as voronoi_sphere_s2.geojson")

from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round,wrap_longitude, wrap_latitude
import geojson

# Initialize RHEALPix DGGS
rdggs = RHEALPixDGGS()

# Specify resolution
resolution = 0
grid = rdggs.grid(resolution)
# print([str(x) for x in grid])
# Function to filter cells crossing the antimeridian
def wrap_latlong(boundary):
    for lat, lon in boundary:
        wrapped_boundary = []
        for lat, lon in boundary:
            wrapped_lat = wrap_latitude(lat)
            wrapped_lon = wrap_longitude(lon)
            wrapped_boundary.append((wrapped_lat, wrapped_lon))
    return wrapped_boundary


# Function to convert cell vertices to GeoJSON Polygon
def cell_to_geojson(cell):
    # Extract vertices and convert to regular floats
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    # print(vertices)
    # Apply the antimeridian filter
    # vertices = wrap_latlong(vertices)
    
    # GeoJSON Polygon requires the last vertex to repeat the first
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    
    # Return GeoJSON feature
    return geojson.Feature(
        geometry=geojson.Polygon([vertices]),
        properties={"rhealpix": str(cell)}  # Add cell ID as a property
    )

# # Convert all cells to GeoJSON features
features = [cell_to_geojson(cell) for cell in grid]
# print(features)
# Create a GeoJSON FeatureCollection
feature_collection = geojson.FeatureCollection(features)

# Save to a GeoJSON file
with open(f"rhealpix_grid_{resolution}.geojson", "w") as f:
    geojson.dump(feature_collection, f, indent=2)

print(f"GeoJSON saved as rhealpix_grid_{resolution}.geojson")

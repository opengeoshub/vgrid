import json
from vgrid.utils import s2
from vgrid.geocode.geocode2geojson import *
from vgrid.utils.s2 import LatLng, CellId, Cell

# latitude, longitude = 10.775275567242561, 106.70679737574993

# s2_resolution = 21 #[0 -->30]
# lat_lng = LatLng.from_degrees(latitude, longitude)
# cell_id = CellId.from_lat_lng(lat_lng)
# cell_id = cell_id.parent(s2_resolution)
# cell_id_token= CellId.to_token(cell_id)
# print(cell_id_token)

# data = s22geojson(cell_id_token)
# print(data)

# output_file = f's2_{s2_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# print(f'GeoJSON written to {output_file}')


import json
from vgrid.utils.s2 import CellId, Cell, LatLng, LatLngRect, RegionCoverer

def s2_covering_to_geojson(min_lat, min_lng, max_lat, max_lng, max_cells):
    # Define the region as a LatLngRect
    region = LatLngRect(
        LatLng.from_degrees(min_lat, min_lng),
        LatLng.from_degrees(max_lat, max_lng),
    )
    
    # Create a RegionCoverer
    coverer = RegionCoverer()
    coverer.min_level = 1
    coverer.max_level = 1
    coverer.max_cells = max_cells
    
    # Get the covering as S2 CellIds
    covering = coverer.get_covering(region)
    
    # Convert S2 Cells to GeoJSON features
    features = []
    for cell_id in covering:
        cell = Cell(cell_id)  # Get the Cell object from CellId
        cell_vertices = []
        
        # Get the vertices of the cell
        for i in range(4):  # Each cell has 4 vertices
            vertex = cell.get_vertex(i)
            latlng = LatLng.from_point(vertex)
            cell_vertices.append([latlng.lng().degrees, latlng.lat().degrees])
        
        # Close the polygon by repeating the first vertex
        cell_vertices.append(cell_vertices[0])
        
        # Add S2 token and Cell ID to properties
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [cell_vertices],
            },
            "properties": {
                "s2_token": cell_id.to_token(),  # Get the S2 Token
                "cell_id": str(cell_id.id()),    # Get the Cell ID as string
            },
        }
        features.append(feature)
    
    # Create the GeoJSON object
    geojson = {
        "type": "FeatureCollection",
        "features": features,
    }
    
    return geojson

# Define the region's bounding box
# min_lat, min_lng, max_lat, max_lng = 8.48, 102.17, 23.53, 109.73  # Example bounding box
min_lat, min_lng, max_lat, max_lng = -90.0, -180.0, 90.0, 180.0

max_cells = 20000  # Adjust based on desired resolution

# Generate the GeoJSON
geojson_data = s2_covering_to_geojson(min_lat, min_lng, max_lat, max_lng, max_cells)

# Save the GeoJSON to a file
output_file = "s2_covering.geojson"
with open(output_file, "w") as f:
    json.dump(geojson_data, f, indent=2)

print(f"GeoJSON saved to {output_file}")

#Reference: 
# https://github.com/aaliddell/s2cell, 
# https://medium.com/@claude.ducharme/selecting-a-geo-representation-81afeaf3bf01
# https://github.com/sidewalklabs/s2sphere
# https://github.com/google/s2geometry/tree/master/src/python
# https://github.com/google/s2geometry
# https://gis.stackexchange.com/questions/293716/creating-shapefile-of-s2-cells-for-given-level
# https://s2sphere.readthedocs.io/en/latest/quickstart.html

import json
from vgrid.geocode.s2 import CellId, Cell, LatLng, LatLngRect, RegionCoverer

def s2_covering_to_geojson(min_lat, min_lng, max_lat, max_lng, max_cells):
    # Define the region as a LatLngRect
    region = LatLngRect(
        LatLng.from_degrees(min_lat, min_lng),
        LatLng.from_degrees(max_lat, max_lng),
    )
    
    # Create a RegionCoverer
    coverer = RegionCoverer()
    coverer.min_level = 1
    coverer.max_level = 2
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
            cell_vertices.append((latlng.lat().degrees, latlng.lng().degrees))
        
        # Close the polygon by repeating the first vertex
        cell_vertices.append(cell_vertices[0])
        
       
        # Add S2 token and Cell ID to properties
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [[(lon, lat) for lat, lon in cell_vertices]],  # Flip lat/lon to lon/lat for GeoJSON
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

# Define the region's bounding box for the whole world
min_lat, min_lng, max_lat, max_lng = -90.0, -180.0, 90.0, 180.0
# min_lat, min_lng, max_lat, max_lng = 8.48, 102.17, 23.53, 109.73  # Example bounding box
max_cells = 2000  # Adjust based on desired resolution

# Generate the GeoJSON
geojson_data = s2_covering_to_geojson(min_lat, min_lng, max_lat, max_lng, max_cells)

# Save the GeoJSON to a file
output_file = "s2_covering.geojson"
with open(output_file, "w") as f:
    json.dump(geojson_data, f, indent=2)

print(f"GeoJSON saved to {output_file}")


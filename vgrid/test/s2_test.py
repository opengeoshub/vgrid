import json
from vgrid.utils import s2
s2_token = '31a06f1'
cell_id = s2.CellId.from_token(s2_token)
cell = s2.Cell(cell_id)
centroid_lat, centroid_lon = s2.CellId.to_lat_lng(cell)
print(centroid_lat, centroid_lon)

# def get_s2_covering(geojson, level, compact=False):
#     """Generate S2 cells covering a GeoJSON geometry at the given level."""
#     geometry = json.loads(geojson)
#         # Extract bounding box of the polygon
#     lats, lngs = zip(*[(lat, lng) for ring in geometry["coordinates"] for lng, lat in ring])
#     min_lat, max_lat = min(lats), max(lats)
#     min_lng, max_lng = min(lngs), max(lngs)
    
#     region = s2.LatLngRect.from_point_pair(
#         s2.LatLng.from_degrees(min_lat, min_lng),
#         s2.LatLng.from_degrees(max_lat, max_lng),
#     )

#     # Generate covering S2 cells
#     coverer = s2.RegionCoverer()
#     coverer.min_level = level
#     coverer.max_level = level
#     covering = coverer.get_covering(region)

#     if compact:
#         # Convert to S2CellUnion and normalize
#         cell_union = s2.CellUnion(covering)
#         # cell_union.normalize()
#         return [cell.id() for cell in cell_union.cell_ids()]  # ✅ Use `.cell_ids()`

#     return [cell.id() for cell in covering]

# # Example GeoJSON Polygon
# geojson_polygon = '{"type":"Polygon","coordinates":[[[-122.5,37.7],[-122.4,37.7],[-122.4,37.8],[-122.5,37.8],[-122.5,37.7]]]}'

# # Generate S2 covering at level 10 with compaction
# s2_cells = get_s2_covering(geojson_polygon, level=10, compact=True)

# print("Compact S2 Cells:", s2_cells)


# import json
# from vgrid.utils import s2
# from vgrid.conversion.dggs2geojson import *
# from vgrid.utils.s2 import LatLng, CellId, Cell

# latitude, longitude = 10.775275567242561, 106.70679737574993

# s2_resolution = 21 #[0..30]
# lat_lng = s2.LatLng.from_degrees(latitude, longitude)
# cell_id = s2.CellId.from_lat_lng(lat_lng)
# cell_id = cell_id.parent(s2_resolution)
# cell_id_token= s2.CellId.to_token(cell_id)
# print(f'S2 Cell Token at resolution {s2_resolution}: {cell_id_token}')
# lat_lng = cell_id.to_lat_lng() 
# print(f'Decode {cell_id_token} to WGS84: {lat_lng}')
# print(f'{cell_id_token} to GeoJSON:\n', s22geojson(cell_id_token))


# s2_resolution = 21 #[0 -->30]
# # latitude = point.y
# # longitude = point.x
# lat_lng = s2.LatLng.from_degrees(latitude, longitude)
# cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
# cell_id = cell_id_max_res.parent(s2_resolution)
# cell = s2.Cell(cell_id)
# cell_token = CellId.to_token(cell.id()) # get Cell ID Token, shorter than cell_id.id()
# print(cell)
# print(cell_token)
# # data = s22geojson(cell_id_token)
# # print(data)

# output_file = f's2_{s2_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# print(f'GeoJSON written to {output_file}')


# import json
# from vgrid.utils.s2 import CellId, Cell, LatLng, LatLngRect, RegionCoverer

# def s2_covering_to_geojson(min_lat, min_lng, max_lat, max_lng, max_cells):
#     # Define the region as a LatLngRect
#     region = LatLngRect(
#         LatLng.from_degrees(min_lat, min_lng),
#         LatLng.from_degrees(max_lat, max_lng),
#     )
    
#     # Create a RegionCoverer
#     coverer = RegionCoverer()
#     coverer.min_level = 1
#     coverer.max_level = 1
#     coverer.max_cells = max_cells
    
#     # Get the covering as S2 CellIds
#     covering = coverer.get_covering(region)
    
#     # Convert S2 Cells to GeoJSON features
#     features = []
#     for cell_id in covering:
#         cell = Cell(cell_id)  # Get the Cell object from CellId
#         cell_vertices = []
        
#         # Get the vertices of the cell
#         for i in range(4):  # Each cell has 4 vertices
#             vertex = cell.get_vertex(i)
#             latlng = LatLng.from_point(vertex)
#             cell_vertices.append([latlng.lng().degrees, latlng.lat().degrees])
        
#         # Close the polygon by repeating the first vertex
#         cell_vertices.append(cell_vertices[0])
        
#         # Add S2 token and Cell ID to properties
#         feature = {
#             "type": "Feature",
#             "geometry": {
#                 "type": "Polygon",
#                 "coordinates": [cell_vertices],
#             },
#             "properties": {
#                 "s2_token": cell_id.to_token(),  # Get the S2 Token
#                 "cell_id": str(cell_id.id()),    # Get the Cell ID as string
#             },
#         }
#         features.append(feature)
    
#     # Create the GeoJSON object
#     geojson = {
#         "type": "FeatureCollection",
#         "features": features,
#     }
    
#     return geojson

# # Define the region's bounding box
# # min_lat, min_lng, max_lat, max_lng = 8.48, 102.17, 23.53, 109.73  # Example bounding box
# min_lat, min_lng, max_lat, max_lng = -90.0, -180.0, 90.0, 180.0

# max_cells = 20000  # Adjust based on desired resolution

# # Generate the GeoJSON
# geojson_data = s2_covering_to_geojson(min_lat, min_lng, max_lat, max_lng, max_cells)

# # Save the GeoJSON to a file
# output_file = "s2_covering.geojson"
# with open(output_file, "w") as f:
#     json.dump(geojson_data, f, indent=2)

# print(f"GeoJSON saved to {output_file}")

import h3, geojson
from vgrid.geocode.geocode2geojson import *

latitude, longitude = 10.775275567242561, 106.70679737574993

h3_precision = 5
h3_code = h3.latlng_to_cell(latitude, longitude, h3_precision)
h3_decode = h3.cell_to_latlng(h3_code)

print(f'latitude, longitude = {latitude},{longitude}')
print(f'H3 code at precision = {h3_precision}: {h3_code}')
print(f'Decode {h3_code} to WGS84 = {h3_decode}')

bbox = h3.cell_to_boundary(h3_code)
print(bbox)

# data = h32geojson(h3_code)
# output_file = f'h3{h3_precision}.geojson'
# with open(output_file, 'w') as f:
#     geojson.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# print(f'GeoJSON written to {output_file}')


# # Define a polygon in GeoJSON-like format (list of lists of coordinates)
# polygon = [
#     [
#         [-122.4089866999972145, 37.783975720493445],
#         [-122.4089866999972145, 37.786111036503398],
#         [-122.4057911999974659, 37.786111036503398],
#         [-122.4057911999974659, 37.783975720493445],
#         [-122.4089866999972145, 37.783975720493445]  # Closing the polygon
#     ]
# ]

# print( h3.is_valid_cell('85283473fffffff'))

# latlng_to_cell = h3.latlng_to_cell(37.3615593, -122.0553238, 5) 
# print( latlng_to_cell)

# latlng = h3.cell_to_latlng('85283473fffffff')
# print(latlng)

# bbox = h3.cell_to_boundary('85283473fffffff')
# print(bbox)

# h = '8928308280fffff'
# grid_disk = h3.grid_disk(h, 1)
# print(f'grid_disk: {grid_disk}')

# maine = [
#         (45.137, -67.137),
#         (44.810, -66.965),
#         (44.325, -68.033),
#         (43.980, -69.060),
#         (43.684, -70.116),
#         (43.090, -70.646),
#         (43.080, -70.751),
#         (43.220, -70.798),
#         (43.368, -70.982),
#         (43.466, -70.944),
#         (45.305, -71.085),
#         (45.460, -70.660),
#         (45.915, -70.305),
#         (46.693, -70.000),
#         (47.448, -69.237),
#         (47.185, -68.905),
#         (47.355, -68.234),
#         (47.066, -67.790),
#         (45.703, -67.791),
#         (45.137, -67.137),
#     ]

# poly = h3.LatLngPoly(maine)
# cells = h3.h3shape_to_cells(poly, 3)
# shape = h3.cells_to_h3shape(cells)

# print(f'Shape to cells: {cells}')
# print(f'Cells to shape: {shape}')

# # Initialize a list to hold the GeoJSON polygons
# geometries = []

# # Loop through each H3 index and convert to GeoJSON format
# for h3_index in cells:
#     # Get the boundary of the H3 hexagon
#     boundary = h3.cell_to_boundary(h3_index)
    
#     # Create the GeoJSON structure for this hexagon
#     polygon_geojson = {
#         "type": "Polygon",
#         "coordinates": [boundary]  # Wrap the boundary in an outer list
#     }
    
#     # Append the polygon to the geometries list
#     geometries.append(polygon_geojson)

# # Create a GeoJSON FeatureCollection to contain all polygons
# geojson_feature_collection = {
#     "type": "FeatureCollection",
#     "features": []
# }

# # Add each polygon as a Feature to the FeatureCollection
# for polygon in geometries:
#     feature = {
#         "type": "Feature",
#         "geometry": polygon,
#         "properties": {
#             "h3_index": h3_index  # Add the H3 index as a property
#         }
#     }
#     geojson_feature_collection["features"].append(feature)

# # Display the resulting GeoJSON FeatureCollection
# output_file = 'h3.geojson'
# with open(output_file, 'w') as f:
#     geojson.dump(geojson_feature_collection, f, indent=2)  # 'indent' makes the JSON output more readable
# print(f'GeoJSON written to {output_file}')

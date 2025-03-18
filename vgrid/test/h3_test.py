
from vgrid.conversion.dggs2geojson import *
latitude, longitude = 10.775275567242561, 106.70679737574993
print(f'Latitude, Longitude: ({latitude}, {longitude})')
h3_resolution = 13 #[0..15]
h3_code = h3.latlng_to_cell(latitude, longitude, h3_resolution)
h3_decode = h3.cell_to_latlng(h3_code)
print(f'H3 code at resolution {h3_resolution}: {h3_code}')
print(f'Decode {h3_code} to WGS84: {h3_decode}')
print(f'{h3_code} to GeoJSON:\n', h32geojson(h3_code))


# import h3, json
# from vgrid.conversion.cell2geojson import *

# latitude, longitude = 10.775275567242561, 106.70679737574993
# # # latitude, longitude =-79.220986356276, 72.57079775696259
# # h3_resolution = 13 #[0-15]
# # h3_code = h3.latlng_to_cell(latitude, longitude, h3_resolution)
# # h3_decode = h3.cell_to_latlng(h3_code)

# # print(f'latitude, longitude = {latitude},{longitude}')
# # print(f'H3 code at resolution = {h3_resolution}: {h3_code}')
# # print(f'Decode {h3_code} to WGS84 = {h3_decode}')
# # data = h32geojson(str(h3_code))
# # print(data)

# # output_file = f'h3_{h3_resolution}.geojson'
# # with open(output_file, 'w') as f:
# #     json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# # print(f'GeoJSON written to {output_file}')
# poly =h3.LatLngPoly(
#          [(37.68, -122.54), (37.68, -122.34), (37.82, -122.34),
#           (37.82, -122.54)],
#      )
# print(h3.cell_to_boundary('862830827ffffff'))
# # print (h3.geo_to_cells(poly, 5))
# # h3_shape  = h3.geo_to_h3shape(poly)
# # geo = "{'type': 'Polygon', 'coordinates': (((-122.54, 37.68), (-122.34, 37.68), (-122.34, 37.82), (-122.54, 37.82), (-122.54, 37.68)),)}"
# # geo_to_h3_shape = h3.geo_to_h3shape(geo)
# # print(geo_to_h3_shape)

# print(h3.grid_path_cells('862830827ffffff','86283095fffffff'))

# # print(h3.cell_to_parent('85283097fffffff'))

# # print(h3.get_base_cell_number('862830827ffffff'))
# # lat_lon_shape = h3.cells_to_h3shape(['862830827ffffff', '862830947ffffff', '862830877ffffff', '862830807ffffff', '862830957ffffff', '86283082fffffff', '86283095fffffff'])
# # geo = h3.h3shape_to_geo(lat_lon_shape)
# # print(geo)
# # print(h3.grid_distance('862830827ffffff', '862830947ffffff'))

# # print(h3.grid_ring('862830827ffffff', k=1))


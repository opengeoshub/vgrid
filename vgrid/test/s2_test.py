import s2, geojson
from vgrid.geocode import s2sphere
from vgrid.geocode.geocode2geojson import *
from vgrid.geocode.s2sphere import LatLng, CellId

latitude, longitude = 10.775275567242561, 106.70679737574993

s2_precision = 21
# lat_lng = s2sphere.LatLng.from_degrees(latitude, longitude)

# Get the S2 CellId for a specific level (e.g., level 10)
# s2_code = s2sphere.CellId.from_lat_lng(lat_lng)
latlng = LatLng.from_degrees(latitude, longitude)

# Get the S2 cell ID at the specified resolution
cell_id = CellId.from_lat_lng(latlng)

s2_code = cell_id.parent(s2_precision)

# print("Cell ID:", cell_id.id())
print(f'latitude, longitude = {latitude},{longitude}')
print(f's2 code at precision = {s2_precision}: {s2_code.id()}')

cell = s2sphere.Cell(s2_code)
center = cell.get_center()
print("Center:", center)

vertices = [cell.get_vertex(i) for i in range(4)]
print (type(vertices))
print (vertices)

data = s22geojson(s2_code)
print(data)
output_file = f's2_{s2_precision}.geojson'
with open(output_file, 'w') as f:
    geojson.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
import h3, json
from vgrid.geocode.geocode2geojson import *

latitude, longitude = 10.775275567242561, 106.70679737574993
# latitude, longitude =-79.220986356276, 72.57079775696259
h3_resolution = 13 #[0-15]
h3_code = h3.latlng_to_cell(latitude, longitude, h3_resolution)
h3_decode = h3.cell_to_latlng(h3_code)

print(f'latitude, longitude = {latitude},{longitude}')
print(f'H3 code at resolution = {h3_resolution}: {h3_code}')
print(f'Decode {h3_code} to WGS84 = {h3_decode}')
data = h32geojson(str(h3_code))
print(data)

output_file = f'h3_{h3_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
from vgrid.geocode import vcode
from vgrid.geocode.geocode2geojson import *
latitude, longitude = 10.775275567242561, 106.70679737574993

vcode_zoom = 23
vcode_code = vcode.latlon2vcode(latitude, longitude, vcode_zoom)
vcode_encode = vcode.vcode2latlon(vcode_code)
print(f'Vcode at zoom level = {vcode_zoom}: {vcode_code}')
print(f'Convert {vcode_code} to WGS84 = {vcode_encode}')

data = vcode.vcode2geojson(vcode_code)
output_file = 'vcode.geojson'
with open(output_file, 'w') as f:
    geojson.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
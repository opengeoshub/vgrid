from vgrid.geocode import vcode
import json
latitude, longitude = 10.775275567242561, 106.70679737574993

resolution = 23
vcode_code = vcode.latlon2vcode(latitude, longitude, resolution)
vcode_encode = vcode.vcode2latlon(vcode_code)
print(f'Vcode at zoom level = {resolution}: {vcode_code}')
print(f'Convert {vcode_code} to WGS84 = {vcode_encode}')

data = vcode.vcode2geojson(vcode_code)
print(data)

output_file = f'vcode_{resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
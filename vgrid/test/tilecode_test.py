from vgrid.conversion.cell2geojson import *
from vgrid.utils import tilecode
import json
latitude, longitude = 10.775275567242561, 106.70679737574993

resolution = 23
tile_code = tilecode.latlon2tilecode(latitude, longitude, resolution)
tile_encode = tilecode.tilecode2latlon(tile_code)
print(f'Tilecode at zoom level = {resolution}: {tile_code}')
print(f'Convert {tile_code} to WGS84 = {tile_encode}')

data = tilecode2geojson(tile_code)
print(data)

output_file = f'tilecode_{resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
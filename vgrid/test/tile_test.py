from vgrid.conversion.cell2geojson import *
from vgrid.utils import tile
import json
latitude, longitude = 10.775275567242561, 106.70679737574993


tile_esolution = 23  # [0..26]
tile_code = tile.latlon2tilecode(latitude, longitude, tile_esolution)
tile_encode = tile.tilecode2latlon(tile_code)
print(f'Tilecode at zoom level = {tile_esolution}: {tile_code}')
print(f'Convert {tile_code} to WGS84 = {tile_encode}')
print(f'{tile_code} to GeoJSON:\n', tilecode2geojson(tile_code))


# resolution = 23
# tile_code = tile.latlon2tilecode(latitude, longitude, resolution)
# tile_encode = tile.tilecode2latlon(tile_code)
# print(f'Tilecode at zoom level = {resolution}: {tile_code}')
# print(f'Convert {tile_code} to WGS84 = {tile_encode}')

# data = tilecode2geojson(tile_code)
# print(data)

# output_file = f'tilecode_{resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# print(f'GeoJSON written to {output_file}')
from vgrid.utils import olc
from vgrid.conversion.dggs2geojson import *
latitude, longitude = 10.775275567242561, 106.70679737574993


olc_resolution = 15 #[10..15]
olc_code = olc.encode(latitude, longitude, olc_resolution)
olc_decode = olc.decode(olc_code)
# print(f'OLC at resolution = {olc_resolution}: {olc_code}')
# print(f'Decode {olc_code} to center and cell in WGS84 = {olc_decode}')
# print(f'{olc_code} to GeoJSON:\n', olc2geojson(olc_code))

coord = olc.decode(olc_code)
msg = (coord.latitudeCenter, coord.longitudeCenter)
print(msg)

# olc_resolution = 15 #[10-->15]
# olc_code = olc.encode(latitude, longitude, olc_resolution)
# # olc_decode = olc.decode(olc_code)
# print(f'OLC at resolution = {olc_resolution}: {olc_code}')
# # print(f'Decode {olc_code} to center and cell in WGS84 = {olc_decode}')

# # data = olc2geojson(olc_code)
# # print(data)
# # output_file = f'olc{olc_resolution}.geojson'
# # with open(output_file, 'w') as f:
# #     json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# # print(f'GeoJSON written to {output_file}')
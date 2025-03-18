from vgrid.utils import mgrs
from vgrid.conversion.dggs2geojson import *
import json
latitude, longitude = 40.00194441, 23.99972080
latitude, longitude = 10.775275567242561, 106.70679737574993

mgrs_resolution = 4 # [0..5]
mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_resolution)
mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
print(f'MGRS Code at resolution {mgrs_resolution}: {mgrs_code}')
print(f'Convert {mgrs_code} to WGS84: {mgrs_code_to_wgs}')
print(f'{mgrs_code} to GeoJSON:\n', mgrs2geojson(mgrs_code))


# mgrs_resolution = 4 # [0 -->5]
# mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_resolution)
# mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
# print(f'MGRS Code at resolution = {mgrs_resolution}: {mgrs_code}')
# print(f'Convert {mgrs_code} to WGS84 = {mgrs_code_to_wgs}')

# data = mgrs2geojson(mgrs_code)
# print(data)
# output_file = f'mgrs{mgrs_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
# print(f'GeoJSON written to {output_file}')
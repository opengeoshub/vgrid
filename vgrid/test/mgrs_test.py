from vgrid.geocode import mgrs
from vgrid.geocode.geocode2geojson import *
latitude, longitude = 40.04042917,18.03399178
latitude, longitude = 40.12268902,18.12483604
latitude, longitude = 10.12521164,102.22014118
latitude, longitude = 10.775275567242561, 106.70679737574993

mgrs_precision = 4
mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_precision)
mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
print(f'MGRS Code at precision = {mgrs_precision}: {mgrs_code}')
print(f'Convert {mgrs_code} to WGS84 = {mgrs_code_to_wgs}')

# mgrs_code_to_wgs_mgrs_code = mgrs.toMgrs(mgrs_code_to_wgs[0], mgrs_code_to_wgs[1], mgrs_precision)
# print(f'Convert {mgrs_code_to_wgs} back to MGRS = {mgrs_code_to_wgs_mgrs_code}')

data = mgrs2geojson(mgrs_code,latitude, longitude)
output_file = f'mgrs{mgrs_precision}.geojson'
with open(output_file, 'w') as f:
    geojson.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')

# mgrs_code = '47PQM'
# # mgrs_code = '48PSS'
# print(mgrs_code)
# zone, hemisphere, easting, northing = mgrs._mgrsToUtm(mgrs_code)
# print(f'MGRS to UTM {zone}, {hemisphere}, {easting}, {northing}')

# print(mgrs.toWgs(mgrs_code))
# print(mgrs._breakMgrsString(mgrs_code))
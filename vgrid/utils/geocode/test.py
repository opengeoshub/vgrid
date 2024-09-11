from vgrid.utils.geocode import mgrs, maidenhead, geohash, olc
latitude, longitude = 10.775275567242561, 106.70679737574993
print(f'latitude, longitude: ({latitude}, {longitude})')

# MGRS
print('\nMGRS:')
mgrs_precision = 5
mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_precision)
mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
print(f'MGRS Code at precision = {mgrs_precision}: {mgrs_code}')
print(f'Convert {mgrs_code} back to WGS84 = {mgrs_code_to_wgs}')

# Maidenhead
print('\nMaidenhead:')
maidenhead_precision = 3
maidenhead_code = maidenhead.toMaiden(latitude, longitude, maidenhead_precision)
maidenGrid = maidenhead.maidenGrid(maidenhead_code)
print(f'Maidenhead Code at precision = {maidenhead_precision}: {maidenhead_code}')
print(f'Convert {maidenhead_code} to center and cell in WGS84 = {maidenGrid}')

# Geohash
print('\nGeohash:')
geohash_precision = 12
geohash_code = geohash.encode(latitude, longitude, geohash_precision)
geohash_decode = geohash.decode(geohash_code, False)
print(f'Geohash Code at precision = {geohash_precision}: {geohash_code}')
print(f'Decode {geohash_code} to center in WGS84 = {geohash_decode}')

# OLC
print('\nOpen Location Code (OLC):')
olc_precision = 10
olc_code = olc.encode(latitude, longitude, olc_precision)
olc_decode = olc.decode(olc_code)
print(f'OLC at precision = {olc_precision}: {olc_code}')
print(f'Decode {olc_code} to center and cell in WGS84 = {olc_decode}')
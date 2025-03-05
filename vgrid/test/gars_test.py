
# from vgrid.utils.gars.garsgrid import GARSGrid
# from vgrid.conversion.cell2geojson import *
# import json 

# latitude, longitude = 10.775275567242561, 106.70679737574993# GARS encoding
# gars_resolution = 1 # 1, 5, 15, 30 minutes
# gars_grid = GARSGrid.from_latlon(latitude, longitude, gars_resolution)
# gars_code = gars_grid.gars_id
# print(gars_code)

# data = gars2geojson(gars_code)
# print(data)
# output_file = f'gars_{gars_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  
# print(f'GeoJSON written to {output_file}')

from vgrid.utils import s2, olc, geohash, georef, mgrs, maidenhead, gars, tilecode
import h3, json
from vgrid.conversion.cell2geojson import *
from vgrid.conversion.latlon2cell import *

latitude, longitude = 10.775275567242561, 106.70679737574993

gars_resolution = 1 # [1, 5, 15, 30 minutes]
gars_grid = gars.garsgrid.GARSGrid.from_latlon(latitude, longitude, gars_resolution)
gars_code = gars_grid.gars_id
print(f'GARS code at resolution {gars_resolution}: {gars_code}')
print(f'{gars_code} to GeoJSON:\n', gars2geojson(gars_code))

# from vgrid.utils import s2, olc, geohash, georef, mgrs, maidenhead, gars
# print('\nGARS:')
# gars_resolution = 1 # [1, 5, 15, 30 minutes]
# gars_grid = gars.garsgrid.GARSGrid.from_latlon(latitude, longitude, gars_resolution)
# gars_code = gars_grid.gars_id
# print(gars_code)

# data = gars2geojson(gars_code)
# output_file = f'gars_{gars_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  
# print(f'GeoJSON written to {output_file}')

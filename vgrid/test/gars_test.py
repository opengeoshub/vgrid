
from vgrid.geocode.gars import GARSGrid
from vgrid.geocode.geocode2geojson import *

latitude, longitude = 10.775275567242561, 106.70679737574993# GARS encoding
# from latlon
gars_precision = 1 # 1, 5, 15, 30 minutes
gars_code = GARSGrid.from_latlon(latitude, longitude, gars_precision)
print(gars_code)

data = gars2geojson(gars_code)
output_file = f'gars_{gars_precision}.geojson'
with open(output_file, 'w') as f:
    geojson.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')

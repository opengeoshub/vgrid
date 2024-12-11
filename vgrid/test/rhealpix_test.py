from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003,UNIT_003
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.geocode.geocode2geojson import *

E = WGS84_ELLIPSOID
# rdggs = RHEALPixDGGS()
# rdggs = WGS84_003
rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
latitude, longitude = 10.775275567242561, 106.70679737574993
resolution = 1
# p = (longitude, latitude)
# rhealpix_resolution = 14 # [0,15]
# rhealpix_code = rdggs.cell_from_point(rhealpix_resolution, p, plane=False)
# print (rhealpix_code)

# geojson_features = rhealpix2geojson(rhealpix_code)
# print(geojson_features)

# output_file = f'rhealpix_{rhealpix_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(geojson_features, f, indent=2)  
# print(f'GeoJSON written to {output_file}')

rhealpix_grid = rhealpix_dggs.grid(resolution)
cell_count = 0

# Iterate over and print each cell
for cell in rhealpix_grid:
    print(cell)
    cell_count += 1  # Increment the count

# Print the total count
print('Number of cells:', cell_count)
print (rhealpix_dggs.num_cells(resolution))
from vgrid.utils.rhealpixdggs.dggs import Cell, RHEALPixDGGS, WGS84_003,UNIT_003,WGS84_003_RADIANS, array
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.conversion.cell2geojson import *

E = WGS84_ELLIPSOID
# rdggs = RHEALPixDGGS()
# rdggs = WGS84_003
rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
latitude, longitude = 10.775275567242561, 106.70679737574993
resolution = 1
# Function to convert degrees to radians

p = (longitude, latitude)
rhealpix_resolution = 2 # [0,15]
rhealpix_cell = rhealpix_dggs.cell_from_point(rhealpix_resolution, p, plane=False)
# rhealpix_cell = rhealpix_dggs.cell(['R',0])
print (rhealpix_cell)
for _, cell in rhealpix_cell.neighbors(plane=False).items():
    print (cell)
    
# rhealpix_cell_west = rhealpix_cell.neighbor('west',plane = False)
# rhealpix_cell_west_north = rhealpix_cell_west.neighbor('north',plane = False)
# rhealpix_cell__west_south = rhealpix_cell_west.neighbor('south',plane = False)
# print(rhealpix_cell_west)
# print(rhealpix_cell_west_north)
# print(rhealpix_cell__west_south)


# rhealpix_cell_east = rhealpix_cell.neighbor('east',plane = False)
# rhealpix_cell_east_north = rhealpix_cell_east.neighbor('north',plane = False)
# rhealpix_cell_east_south = rhealpix_cell_east.neighbor('south',plane = False)
# print(rhealpix_cell_east)
# print(rhealpix_cell_east_north)
# print(rhealpix_cell_east_south)


# rhealpix_cell_north = rhealpix_cell.neighbor('north',plane = False)
# print(rhealpix_cell_north)

# rhealpix_cell_south = rhealpix_cell.neighbor('south',plane = False)
# print(rhealpix_cell_south)


# print(rhealpix_cell.centroid(plane = False))
# print(rhealpix_cell.nucleus(plane = False))
# print(rhealpix_cell.color.__name__)
# parent_cell = rhealpix_dggs.cell_from_point(rhealpix_resolution-1, p, plane=False)
# print (parent_cell)
# print([str(cell) for cell in parent_cell.subcells()])

# for (direction, cell) in sorted(rhealpix_cell.neighbors(plane=True).items()):
#    print(direction, cell)

# print(rhealpix_dggs.triangle(*rhealpix_cell.nucleus(), inverse=True))


    
# print(rhealpix_cell.neighbor('west',plane=False))

    
# geojson_features = rhealpix2geojson(rhealpix_code)
# print(geojson_features)

# output_file = f'rhealpix_{rhealpix_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(geojson_features, f, indent=2)  
# print(f'GeoJSON written to {output_file}')

# rhealpix_grid = rhealpix_dggs.grid(resolution)
# cell_count = 0

# # Iterate over and print each cell
# for cell in rhealpix_grid:
#     print(cell)
#     cell_count += 1  # Increment the count

# # Print the total count
# print('Number of cells:', cell_count)
# print (rhealpix_dggs.num_cells(resolution))
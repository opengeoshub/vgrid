from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003,UNIT_003,WGS84_003_RADIANS, array
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.geocode.geocode2geojson import *

E = WGS84_ELLIPSOID
# rdggs = RHEALPixDGGS()
# rdggs = WGS84_003
rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
latitude, longitude = 10.775275567242561, 106.70679737574993
resolution = 1
# Function to convert degrees to radians


def deg_to_rad(degrees):
    return degrees * (math.pi / 180)

# Define the bounding box coordinates in degrees
ul = (106.6908, 10.787395)  # Upper-left corner (longitude, latitude)
dr = (106.733069, 10.759213)  # Lower-right corner (longitude, latitude)
rhealpix_dggs = UNIT_003

rdggs = WGS84_003_RADIANS
R_A = rdggs.ellipsoid.R_A
# Convert the coordinates to radians (if necessary)
ul_rad = (deg_to_rad(ul[0]), deg_to_rad(ul[1]))
dr_rad = (deg_to_rad(dr[0]), deg_to_rad(dr[1]))
pi = 3.14
p = (0, pi/12)
q = (pi/6 - 1e-6, 0)
ul = R_A*array((-0.1, pi/4))
dr = R_A*array((0.1, -pi/4))  # Rectangle
            

 
# Assuming `dggs` is your DGGS grid system object
result = rdggs.cells_from_region(1, p, q , plane=False)  # Assuming you are working with planar DGGS
for row in result:
    print([str(cell) for cell in row])
 
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

# rhealpix_grid = rhealpix_dggs.grid(resolution)
# cell_count = 0

# # Iterate over and print each cell
# for cell in rhealpix_grid:
#     print(cell)
#     cell_count += 1  # Increment the count

# # Print the total count
# print('Number of cells:', cell_count)
# print (rhealpix_dggs.num_cells(resolution))
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003,UNIT_003
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.utils.rhealpixdggs.utils import my_round

E = WGS84_ELLIPSOID
# rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=2, N_side=3)
# rdggs = RHEALPixDGGS()
rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
# rdggs = WGS84_003
# print(rdggs)
latitude, longitude = 10.775275567242561, 106.70679737574993
p = (longitude, latitude)
resolution = 0 # [0,15]
c = rdggs.cell_from_point(15, p, plane=False)
# rdggs.cell_area(1)
# print(rdggs.cell_area(1))
print (c)
for p in c.vertices(plane=False, trim_dart=True):
    print(tuple(x.tolist() for x in p))

# c = rdggs.cell(['N', 7])
# print(rdggs.triangle(*c.nucleus(), inverse=False))

# c = rdggs.cell(['R', 8])
# cell = rdggs.cell(('R3126036255355',))
# for p in c.vertices():
# for p in c.vertices(plane=False, trim_dart=True):
#     print(tuple(x.tolist() for x in p))
# for p in c.vertices():
#     print (p)
# print(c)
# def rhealpix_cell_from_id(cell_id):
#     # Example parsing logic (adjust according to your format)
#     region = int(cell_id[1:2])  # Extract the region (e.g., 'R3' -> 3)
#     level = int(cell_id[3:4])   # Extract the level (e.g., 'L2' -> 2)
#     pixel = int(cell_id[5:])    # Extract the pixel index (e.g., 'P12345' -> 12345)

#     # Return a dictionary or a custom object representing the cell
#     return {
#         "region": region,
#         "level": level,
#         "pixel": pixel
#     }

# Example usage
# cell_id = "N"
# cell = rhealpix_cell_from_id(cell_id)
# def rhealpix_to_rdggs(rhealpix_cell):
#     """
#     Convert a RHEALPix cell to an RDGGS cell.

#     Args:
#         rhealpix_cell (dict): Dictionary with keys 'region', 'level', 'pixel'.
#                               Example: {'region': 3, 'level': 2, 'pixel': 12345}

#     Returns:
#         rdggs.Cell: An RDGGS cell object.
#     """
#     # Map RHEALPix components to RDGGS components
#     region = rhealpix_cell['region']
#     resolution = rhealpix_cell['level']  # Assuming 'level' corresponds to RDGGS resolution
#     pixel_index = rhealpix_cell['pixel']
#     print('region:', region)
#     # Create an RDGGS cell
#     rdggs_cell = rdggs.cell((region, resolution, pixel_index))

#     return rdggs_cell

# cell = rhealpix_to_rdggs(cell)
# print(cell)
# for p in cell.vertices():
#     print (p)


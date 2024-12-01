from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003,UNIT_003
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.utils.rhealpixdggs.utils import my_round

# E = WGS84_ELLIPSOID
# rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=2, N_side=3)
# rdggs = RHEALPixDGGS(N_side=15)
rdggs = RHEALPixDGGS()
# rdggs = WGS84_003
print(rdggs)
latitude, longitude = 10.775275567242561, 106.70679737574993
p = (longitude, latitude)
c = rdggs.cell_from_point(13, p, plane=False)
# c = rdggs.cell(['R', 8])
print(c)
# for p in c.vertices():
for p in c.vertices(plane=False, trim_dart=True):
    print(tuple(x.tolist() for x in p))
# for p in c.vertices():
#     print (p)

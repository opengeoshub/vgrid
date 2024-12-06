from vgrid.geocode import mgrs, maidenhead, geohash, georef, olc, s2, tilecode
import h3
from vgrid.geocode.s2 import LatLng, CellId

from vgrid.utils.gars.garsgrid import GARSGrid
from vgrid.utils import mercantile

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.utils import my_round
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.utils.eaggr.enums.model import Model

import json, re
import argparse

def latlon2h3(lat,lon,res=13):
    # res: [0..15]
    h3_cell = h3.latlng_to_cell(lat, lon, res)
    return h3_cell

def latlon2h3_cli():
    """
    Command-line interface for latlon2h3.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to H3 code at a specific Resolution [0.15]. \
                                     Usage: latlon2h3 <lat> <lon> <res>. \
                                     Ex: latlon2h3 10.775275567242561, 106.70679737574993, 13")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..15]")
    args = parser.parse_args()
    h3_cell = latlon2h3(args.lat,args.lon,args.res)
    print(h3_cell)

def latlon2s2(lat,lon,res=21):
    # res: [0..30]
    lat_lng = LatLng.from_degrees(lat, lon)
    cell_id = CellId.from_lat_lng(lat_lng)
    cell_id = cell_id.parent(res)
    cell_token = CellId.to_token(cell_id)
    return cell_token

def latlon2s2_cli():
    """
    Command-line interface for latlon2s2.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to S2 code at a specific Resolution [0..30]. \
                                     Usage: latlon2s2 <lat> <lon> <res>. \
                                     Ex: latlon2s2 10.775275567242561 106.70679737574993 21")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..30]")
    args = parser.parse_args()
    s2_cell = latlon2s2(args.lat,args.lon,args.res)
    print(s2_cell)

def latlon2rhealpix(lat,lon,res=14):
    # res: [0..15]
    E = WGS84_ELLIPSOID
    rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
    point = (lon, lat)
    rhealpix_cell = rdggs.cell_from_point(res, point, plane=False)
    return rhealpix_cell

def latlon2rhealpix_cli():
    """
    Command-line interface for latlon2rhealpix.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Rhealpix code at a specific Resolution [0..15]. \
                                     Usage: latlon2rhealpix <lat> <lon> <res>. \
                                     Ex: latlon2rhealpix 10.775275567242561 106.70679737574993 14")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..15]")
    args = parser.parse_args()
    rhealpix_cell = latlon2rhealpix(args.lat,args.lon,args.res)
    print(rhealpix_cell)

def latlon2eaggrisea4t(lat,lon,res=21):
    # res: [0..38]
    eaggr_dggs = Eaggr(Model.ISEA4T)
    max_accuracy = 10**(-10) # maximum cell_id length with 40 characters
    lat_long_point = LatLongPoint(lat, lon, max_accuracy)
    eaggr_cell_max_accuracy = eaggr_dggs.convert_point_to_dggs_cell(lat_long_point)
    cell_id_len = res+2
    eaggr_cell = DggsCell(eaggr_cell_max_accuracy._cell_id[:cell_id_len])
    return eaggr_cell._cell_id

def latlon2eaggrisea4t_cli():
    """
    Command-line interface for latlon2eaggrisea4t.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to EaggrISEA4T code at a specific Resolution [0..38]. \
                                     Usage: latlon2eaggrisea4t <lat> <lon> <res>. \
                                     Ex: latlon2eaggrisea4t 10.775275567242561 106.70679737574993 21")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..38]")
    args = parser.parse_args()
    eaggr_cell = latlon2eaggrisea4t(args.lat,args.lon,args.res)
    print(eaggr_cell)


def latlon2olc(lat,lon,res=11):
    # res: [10..15]
    olc_cell = olc.encode(lat, lon, res)
    return olc_cell

def latlon2olc_cli():
    """
    Command-line interface for latlon2olc.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to OLC/ Google Plus Code at a specific Code length [10..15]. \
                                     Usage: latlon2olc <lat> <lon> <res>. \
                                     Ex: latlon2olc 10.775275567242561 106.70679737574993 11")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution/ Code length [10..15]")
    args = parser.parse_args()
    olc_cell = latlon2olc(args.lat,args.lon,args.res)
    print(olc_cell)

def latlon2geohash(lat,lon,res=9):
    # res: [1..30]
    geohash_cell = geohash.encode(lat, lon, res)
    return geohash_cell

def latlon2geohash_cli():
    """
    Command-line interface for latlon2geohash.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Geohash code at a specific resolution [1..30]. \
                                     Usage: latlon2geohash <lat> <lon> <res>. \
                                     Ex: latlon2geohash 10.775275567242561 106.70679737574993 9")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1..30]")
    args = parser.parse_args()
    geohash_cell = latlon2geohash(args.lat,args.lon,args.res)
    print(geohash_cell)

def latlon2georef(lat,lon,res=4):
    # res: [0..10]
    georef_cell = georef.encode(lat,lon,res)
    return georef_cell

def latlon2georef_cli():
    """
    Command-line interface for latlon2georef.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to GEOREF code at a specific resolution [0..10]. \
                                     Usage: latlon2georef <lat> <lon> <res>. \
                                     Ex: latlon2georef 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..10]")
    args = parser.parse_args()
    georef_cell = latlon2georef(args.lat,args.lon,args.res)
    print(georef_cell)

def latlon2mgrs(lat,lon,res=4):
    # res: [0..5]
    mgrs_cell = mgrs.toMgrs(lat,lon,res)
    return mgrs_cell

def latlon2mgrs_cli():
    """
    Command-line interface for latlon2mgrs.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to GEOREF code at a specific resolution [0..5]. \
                                     Usage: latlon2mgrs <lat> <lon> <res>. \
                                     Ex: latlon2mgrs 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution  [0..5]")
    args = parser.parse_args()
    mgrs_cell = latlon2mgrs(args.lat,args.lon,args.res)
    print(mgrs_cell)

def latlon2tilecode(lat,lon,res=23):
    # res: [0..26]
    tilecode_cell = tilecode.latlon2tilecode(lat,lon,res)
    return tilecode_cell

def latlon2tilecode_cli():
    """
    Command-line interface for latlon2tilecode.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Tile code at a specific resolution/ zoom level [0..26]. \
                                     Usage: latlon2tilecode <lat> <lon> <res>. \
                                     Ex: latlon2tilecode 10.775275567242561 106.70679737574993 23")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution/ Zoom level [0..26]")
    args = parser.parse_args()
    tilecode_cell = latlon2tilecode(args.lat,args.lon,args.res)
    print(tilecode_cell)


def latlon2maidenhead(lat,lon,res=4):
    # res: [1..4]
    maidenhead_cell = maidenhead.toMaiden(lat,lon,res)
    return maidenhead_cell

def latlon2maidenhead_cli():
    """
    Command-line interface for latlon2maidenhead.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Tile code at a specific resolution [1..4]. \
                                     Usage: latlon2maidenhead <lat> <lon> <res>. \
                                     Ex: latlon2maidenhead 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1..4]")
    args = parser.parse_args()
    maidenhead_cell = latlon2maidenhead(args.lat,args.lon,args.res)
    print(maidenhead_cell)

def latlon2gars(lat,lon,res=1):
    # res: [1, 5, 15, 30 minutes]
    gars_cell = gars_grid = GARSGrid.from_latlon(lat,lon,res)
    return gars_cell

def latlon2gars_cli():
    """
    Command-line interface for latlon2gars.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Tile code at a specific resolution [1, 5, 15, 30 minutes]. \
                                     Usage: latlon2gars <lat> <lon> <res>. \
                                     Ex: latlon2gars 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1, 5, 15, 30 minutes]")
    args = parser.parse_args()
    gars_cell = latlon2gars(args.lat,args.lon,args.res)
    print(gars_cell)
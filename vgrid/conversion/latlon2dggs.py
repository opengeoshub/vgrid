from vgrid.utils import s2, olc, geohash, georef, mgrs, maidenhead, tilecode, qtm
import h3

from vgrid.utils.gars.garsgrid import GARSGrid

from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
import platform, math 

if (platform.system() == 'Windows'):
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.generator.isea3hgrid import isea3h_res_accuracy_dict
    from vgrid.generator.isea4tgrid import isea4t_res_accuracy_dict
    
if (platform.system() == 'Linux'):
    from vgrid.utils.dggrid4py import DGGRIDv7, dggs_types
    import geopandas as gpd
    from vgrid.utils.dggrid4py.dggrid_runner import output_address_types
    
from shapely import Point

from vgrid.utils.easedggs.dggs.grid_addressing import geos_to_grid_ids

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
                                     Usage: latlon2h3 <lat> <lon> <res> [0..15]. \
                                     Ex: latlon2h3 10.775275567242561 106.70679737574993 13")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..15]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 15:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..15].")
        return  
    
    h3_cell = latlon2h3(args.lat,args.lon,args.res)
    print(h3_cell)

def latlon2s2(lat,lon,res=21):
    # res: [0..30] 
    lat_lng = s2.LatLng.from_degrees(lat, lon)
    cell_id = s2.CellId.from_lat_lng(lat_lng) # return S2 cell at max level 30
    cell_id = cell_id.parent(res) # get S2 cell at resolution
    cell_token = s2.CellId.to_token(cell_id) # get Cell ID Token, shorter than cell_id.id()
    return cell_token

def latlon2s2_cli():
    """
    Command-line interface for latlon2s2.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to S2 code at a specific Resolution [0..30]. \
                                     Usage: latlon2s2 <lat> <lon> <res> [0..30]. \
                                     Ex: latlon2s2 10.775275567242561 106.70679737574993 21")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..30]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 30:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..30].")
        return  
    
    s2_cell = latlon2s2(args.lat,args.lon,res)
    print(s2_cell)

def latlon2rhealpix(lat,lon,res=14):
    # res: [0..15]        
    E = WGS84_ELLIPSOID
    rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
    point = (lon, lat)
    rhealpix_cell = rdggs.cell_from_point(res, point, plane=False)
    return str(rhealpix_cell)

def latlon2rhealpix_cli():
    """
    Command-line interface for latlon2rhealpix.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Rhealpix code at a specific Resolution [0..15]. \
                                     Usage: latlon2rhealpix <lat> <lon> <res> [0..15]. \
                                     Ex: latlon2rhealpix 10.775275567242561 106.70679737574993 14")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..15]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 15:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..15].")
        return  
    
    rhealpix_cell = latlon2rhealpix(args.lat,args.lon,res)
    print(rhealpix_cell)

def latlon2isea4t(lat,lon,res=21):
    if (platform.system() == 'Windows'):
        # res: [0..39]
        isea4t_dggs = Eaggr(Model.ISEA4T)
        max_accuracy =  isea4t_res_accuracy_dict[39] # maximum cell_id length with 41 characters
        lat_long_point = LatLongPoint(lat, lon, max_accuracy)
        isea4t_cell_max_accuracy = isea4t_dggs.convert_point_to_dggs_cell(lat_long_point)
        cell_id_len = res+2
        isea4t_cell = DggsCell(isea4t_cell_max_accuracy._cell_id[:cell_id_len])
        return isea4t_cell._cell_id

def latlon2isea4t_cli():
    """
    Command-line interface for latlon2isea4t.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to OpenEaggr ISEA4T code at a specific Resolution [0..39]. \
                                     Usage: latlon2isea4t <lat> <lon> <res> [0..39]. \
                                     Ex: latlon2isea4t 10.775275567242561 106.70679737574993 21")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..39]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 39:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..39].")
        return  
    
    eaggr_cell = latlon2isea4t(args.lat,args.lon,res)
    print(eaggr_cell)


def latlon2isea3h(lat,lon,res=27):
    if (platform.system() == 'Windows'):
        # res: [0..40], res=27 is suitable for geocoding
        isea3h_dggs = Eaggr(Model.ISEA3H)
        accuracy = isea3h_res_accuracy_dict.get(res)            
        lat_long_point = LatLongPoint(lat, lon, accuracy)
        isea3h_cell = isea3h_dggs.convert_point_to_dggs_cell(lat_long_point)
        return isea3h_cell.get_cell_id()

def latlon2isea3h_cli():
    """
    Command-line interface for latlon2isea3h.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to OpenEaggr ISEA3H code at a specific Resolution [0..40]. \
                                     Usage: latlon2isea3h <lat> <lon> <res> [0..40]. \
                                     Ex: latlon2isea3h 10.775275567242561 106.70679737574993 14")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..40]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 40:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..40].")
        return  
        
    isea3h_cell = latlon2isea3h(args.lat,args.lon,res)
    print(isea3h_cell)


def latlon2dggrid(lat,lon,dggs_type,res,address_type='SEQNUM'):
    if (platform.system() == 'Linux'):
        dggrid_instance = DGGRIDv7(executable='/usr/local/bin/dggrid', working_dir='.', capture_logs=False, silent=True, tmp_geo_out_legacy=False, debug=False)
        point = Point(lon, lat)
        geodf_points_wgs84 = gpd.GeoDataFrame([{'geometry': point}], crs="EPSG:4326")
        dggrid_cell =  dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only = True, dggs_type = dggs_type,resolution = res)    
        seqnum = dggrid_cell.loc[0, 'seqnums']
        address_type_transform = dggrid_instance.address_transform([seqnum], dggs_type= dggs_type, resolution = res,mixed_aperture_level=None, input_address_type='SEQNUM', output_address_type=address_type)
        # cell_id_api = {f'{dggs_type}_{res}_{address_type_transform.columns[1]}': address_type_transform.loc[0,address_type]}  
        dggrid_cell_id = address_type_transform.loc[0,address_type]
        return dggrid_cell_id
        # return address_type_transform

def latlon2dggrid_cli():
    """
    Command-line interface for latlon2dggrid.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to DGGRID cell at a specific Resolution. \
                                     Usage: latlon2dggrid <lat> <lon> <dggs_type> <res>. \
                                     Ex: latlon2dggrid  10.775275567242561 106.70679737574993 ISEA7H 13")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument('dggs_type', choices=dggs_types, help="Select a DGGS type from the available options.")
    parser.add_argument("res",type=int, help="Resolution")
    parser.add_argument('address_type', choices=output_address_types, 
                        default='SEQNUM',
                        nargs='?',  # This makes the argument optional
                        help="Select an output address type from the available options.")

    args = parser.parse_args()
    dggs_type = args.dggs_type
    res = args.res
    address_type  = args.address_type
    dggrid_cell_id = latlon2dggrid(args.lat,args.lon,dggs_type,res,address_type)
    print(dggrid_cell_id)


def latlon2ease(lat,lon,res=6):
    # res = [0..6]  
    easedggs_cell = geos_to_grid_ids([(lon,lat)],level = res)
    easedggs_cell_id = easedggs_cell['result']['data'][0]
    return easedggs_cell_id

def latlon2ease_cli():
    """
    Command-line interface for latlon2isea3h.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to EASE-DGGS cell at a specific Resolution [0..6]. \
                                     Usage: latlon2ease <lat> <lon> <res> [0..6]. \
                                     Ex: latlon2ease 10.775275567242561 106.70679737574993 6")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [0..6]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 6:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..6].")
        return  
    
    
    easedggs_cell = latlon2ease(args.lat,args.lon,res)
    print(easedggs_cell)


def latlon2qtm(lat, lon, res=10):
    return qtm.latlon_to_qtm_id(lat, lon, res)
    
def latlon2qtm_cli():
    """
    Command-line interface for latlon2qtm.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to QTM. \
                                     Usage: latlon2qtm <lat> <lon> <res> [1..24]. \
                                     Ex: latlon2qtm 10.775275567242561 106.70679737574993 10")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1..24]")
    args = parser.parse_args()
    
    res = args.res
    if res < 1 or res > 24:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [1..24].")
        return  
    
    qtm_id = latlon2qtm(args.lat,args.lon,res)
    print(qtm_id)
    
 
def latlon2olc(lat,lon,res=11):
    # res: [2,4,6,8,10..15]        
    olc_cell = olc.encode(lat, lon, res)
    return olc_cell

def latlon2olc_cli():
    """
    Command-line interface for latlon2olc.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to OLC/ Google Plus Code at a specific Code length [10..15]. \
                                     Usage: latlon2olc <lat> <lon> <res> [2,4,6,8,10..15]. \
                                     Ex: latlon2olc 10.775275567242561 106.70679737574993 11")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution/ Code length [2,4,6,8,10..15]")
    args = parser.parse_args()
    
    res = args.res
    if res not in [2,4,6,8,10,11,12,13,14,15]:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [2,4,6,8,10..15].")
        return  
    
    olc_cell = latlon2olc(args.lat,args.lon,res)
    print(olc_cell)

def latlon2geohash(lat,lon,res=6):
    # res: [1..10]    
    geohash_id= geohash.encode(lat, lon, res)
    return geohash_id

def latlon2geohash_cli():
    """
    Command-line interface for latlon2geohash.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Geohash code at a specific resolution [1..10]. \
                                     Usage: latlon2geohash <lat> <lon> <res>[1..10]. \
                                     Ex: latlon2geohash 10.775275567242561 106.70679737574993 6")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1..10]")
    args = parser.parse_args()
    
    res = args.res
    if res < 1 or res > 10:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [1..10].")
        return  
    
    geohash_id = latlon2geohash(args.lat,args.lon,res)
    print(geohash_id)

def latlon2georef(lat,lon,res=4):
    # res: [-1..5]        
    georef_cell = georef.encode(lat,lon,res)
    return georef_cell

def latlon2georef_cli():
    """
    Command-line interface for latlon2georef.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to GEOREF code at a specific resolution [0..5]. \
                                     Usage: latlon2georef <lat> <lon> <res> [-1..5]. \
                                     Ex: latlon2georef 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [-1..5]")
    args = parser.parse_args()
    
    res = args.res
    if res < -1 or res > 5:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [-1..5].")
        return  
    
    georef_cell = latlon2georef(args.lat,args.lon,res)
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
                                     Usage: latlon2mgrs <lat> <lon> <res> [0..5]. \
                                     Ex: latlon2mgrs 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution  [0..5]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 5:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..5].")
        return  

    mgrs_cell = latlon2mgrs(args.lat,args.lon,res)
    print(mgrs_cell)

def latlon2tilecode(lat,lon,res=23):
    # res: [0..29]        
    tilecode_id = tilecode.latlon2tilecode(lat,lon,res)
    return tilecode_id

def latlon2tilecode_cli():
    """
    Command-line interface for latlon2tilecode.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Tile code at a specific resolution/ zoom level [0..29]. \
                                     Usage: latlon2tilecode <lat> <lon> <res> [0..29]. \
                                     Ex: latlon2tilecode 10.775275567242561 106.70679737574993 23")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution/ Zoom level [0..29]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 29:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..29].")
        return 
    
    tilecode_cell = latlon2tilecode(args.lat,args.lon,res)
    print(tilecode_cell)


def latlon2quadkey(lat,lon,res=23):
    # res: [0..29]        
    quadkey = tilecode.latlon2quadkey(lat,lon,res)
    return quadkey

def latlon2quadkey_cli():
    """
    Command-line interface for latlon2tilecode.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Quadkey at a specific resolution/ zoom level [0..29]. \
                                     Usage: latlon2quadkey <lat> <lon> <res> [0..29]. \
                                     Ex: latlon2quadkey 10.775275567242561 106.70679737574993 23")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution/ Zoom level [0..29]")
    args = parser.parse_args()
    
    res = args.res
    if res < 0 or res > 29:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in [0..29].")
        return 
    
    quadkey = latlon2quadkey(args.lat,args.lon,res)
    print(quadkey)


def latlon2maidenhead(lat,lon,res=4):  
    maidenhead_cell = maidenhead.toMaiden(lat,lon,res)
    return maidenhead_cell

def latlon2maidenhead_cli():
    """
    Command-line interface for latlon2maidenhead.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Tile code at a specific resolution [1..4]. \
                                     Usage: latlon2maidenhead <lat> <lon> <res> [1..4]. \
                                     Ex: latlon2maidenhead 10.775275567242561 106.70679737574993 4")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1..4]")
    args = parser.parse_args()

    res = args.res
    if res < 1 or res > 4:
        print(f"Error: Invalid resolution {args.res}. Please input a valid resolutions in [1..4].")
        return 
        
    maidenhead_cell = latlon2maidenhead(args.lat,args.lon,res)
    print(maidenhead_cell)

def latlon2gars(lat,lon,res=1):
    # res: [1, 5, 15, 30 minutes]
    gars_cell = GARSGrid.from_latlon(lat,lon,res)
    return gars_cell

def latlon2gars_cli():
    """
    Command-line interface for latlon2gars.
    """
    parser = argparse.ArgumentParser(description="Convert Lat, Long to Tile code at a specific resolution [1, 5, 15, 30 minutes]. \
                                     Usage: latlon2gars <lat> <lon> <res> [1, 5, 15, 30 minutes]. \
                                     Ex: latlon2gars 10.775275567242561 106.70679737574993 1")
    parser.add_argument("lat",type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument("res",type=int, help="Input Resolution [1, 5, 15, 30 minutes]")
    args = parser.parse_args()
    
    res = args.res
    valid_resolutions = [1, 5, 15, 30]
    if args.res not in valid_resolutions:
        print(f"Error: Invalid resolution {res}. Please input a valid resolutions in {valid_resolutions} (minutes).")
        return 
            
    gars_cell = latlon2gars(args.lat,args.lon,res)
    print(gars_cell)
from vgrid.utils import s2, olc, geohash, georef, mgrs, mercantile, maidenhead
from vgrid.utils.gars import garsgrid
from vgrid.utils.qtm import constructGeometry, qtm_id_to_facet
import h3

from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
import platform

if (platform.system() == 'Windows'):   
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.generator.isea4tgrid import fix_isea4t_wkt, fix_isea4t_antimeridian_cells
    from vgrid.generator.isea3hgrid import isea3h_cell_to_polygon

if (platform.system() == 'Linux'):
    from vgrid.utils.dggrid4py import DGGRIDv7, dggs_types
    from vgrid.utils.dggrid4py.dggrid_runner import input_address_types


from vgrid.utils.easedggs.constants import levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos

from shapely.wkt import loads
from shapely.geometry import shape, Polygon,mapping

import json, re,os,argparse
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.rhealpixgrid import fix_rhealpix_antimeridian_cells

from vgrid.utils.antimeridian import fix_polygon

from vgrid.generator.settings import graticule_dggs_to_feature, geodesic_dggs_to_feature,isea3h_accuracy_res_dict

from pyproj import Geod
geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID

def h32geojson(h3_ids):
    """
    Convert a list of H3 cell IDs to a GeoJSON FeatureCollection.
    Accepts a single h3_id (string) or a list of h3_ids.
    Skips invalid or error-prone cells.
    """
    if isinstance(h3_ids, str):
        h3_ids = [h3_ids]
    h3_features = []
    for h3_id in h3_ids:
        try:
            cell_boundary = h3.cell_to_boundary(h3_id)
            if cell_boundary:
                filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
                reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
                cell_polygon = Polygon(reversed_boundary)
                resolution = h3.get_resolution(h3_id)
                num_edges = 6
                if h3.is_pentagon(h3_id):
                    num_edges = 5
                h3_feature = geodesic_dggs_to_feature("h3", h3_id, resolution, cell_polygon, num_edges)
                h3_features.append(h3_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": h3_features
    }


def h32geojson_cli():
    """
    Command-line interface for h32geojson supporting multiple H3 cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert H3 cell ID(s) to GeoJSON")
    parser.add_argument("h3", nargs='+', help="Input H3 cell ID(s), e.g., h32geojson 8d65b56628e46bf 8d65b56628e46c3")
    args = parser.parse_args()
    geojson_data = json.dumps(h32geojson(args.h3))
    print(geojson_data)
    

def s22geojson(s2_tokens):
    """
    Convert a list of S2 cell tokens to a GeoJSON FeatureCollection.
    Accepts a single s2_token (string) or a list of s2_tokens.
    Skips invalid or error-prone cells.
    """
    if isinstance(s2_tokens, str):
        s2_tokens = [s2_tokens]
    s2_features = []
    for s2_token in s2_tokens:
        try:
            cell_id = s2.CellId.from_token(s2_token)
            cell = s2.Cell(cell_id)
            if cell:
                vertices = [cell.get_vertex(i) for i in range(4)]
                shapely_vertices = []
                for vertex in vertices:
                    lat_lng = s2.LatLng.from_point(vertex)
                    longitude = lat_lng.lng().degrees
                    latitude = lat_lng.lat().degrees
                    shapely_vertices.append((longitude, latitude))
                shapely_vertices.append(shapely_vertices[0])
                cell_polygon = fix_polygon(Polygon(shapely_vertices))
                resolution = cell_id.level()
                num_edges = 4
                s2_feature = geodesic_dggs_to_feature("s2", s2_token, resolution, cell_polygon, num_edges)
                s2_features.append(s2_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": s2_features
    }

def s22geojson_cli():
    """
    Command-line interface for s22geojson supporting multiple S2 cell tokens.
    """
    parser = argparse.ArgumentParser(description="Convert S2 cell token(s) to GeoJSON")
    parser.add_argument("s2", nargs='+', help="Input S2 cell token(s), e.g., s22geojson 31752f45cc94 31752f45cc95")
    args = parser.parse_args()
    geojson_data = json.dumps(s22geojson(args.s2))
    print(geojson_data)


def rhealpix_cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = fix_rhealpix_antimeridian_cells(vertices)
    return Polygon(vertices)

def rhealpix2geojson(rhealpix_ids):
    if isinstance(rhealpix_ids, str):
        rhealpix_ids = [rhealpix_ids]
    rhealpix_features = []
    for rhealpix_id in rhealpix_ids:
        try:
            rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
            rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            if rhealpix_cell:
                resolution = rhealpix_cell.resolution
                cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
                num_edges = 4
                if rhealpix_cell.ellipsoidal_shape() == 'dart':
                    num_edges = 3
                rhealpix_feature = geodesic_dggs_to_feature("rhealpix", rhealpix_id, resolution, cell_polygon, num_edges)
                rhealpix_features.append(rhealpix_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": rhealpix_features
    }

def rhealpix2geojson_cli():
    """
    Command-line interface for rhealpix2geojson supporting multiple Rhealpix cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert Rhealpix cell ID(s) to GeoJSON")
    parser.add_argument("rhealpix", nargs='+', help="Input Rhealpix cell ID(s), e.g., rhealpix2geojson R31260335553825 R31260335553826")
    args = parser.parse_args()
    geojson_data = json.dumps(rhealpix2geojson(args.rhealpix))
    print(geojson_data)
    

def isea4t2geojson(isea4t_ids):
    if isinstance(isea4t_ids, str):
        isea4t_ids = [isea4t_ids]
    isea4t_features = []
    if (platform.system() == 'Windows'):
        for isea4t_id in isea4t_ids:
            try:
                isea4t_dggs = Eaggr(Model.ISEA4T)
                cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(isea4t_id), ShapeStringFormat.WKT)
                cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
                if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14') \
                    or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
                    cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
                if cell_to_shape_fixed:
                    resolution = len(isea4t_id) - 2
                    num_edges = 3
                    cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
                    isea4t_feature = geodesic_dggs_to_feature("isea4t", isea4t_id, resolution, cell_polygon, num_edges)
                    isea4t_features.append(isea4t_feature)
            except Exception:
                continue
        return {
            "type": "FeatureCollection",
            "features": isea4t_features
        }

def isea4t2geojson_cli():
    """
    Command-line interface for isea4t2geojson supporting multiple ISEA4T cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert Open-Eaggr ISEA4T cell ID(s) to GeoJSON")
    parser.add_argument("isea4t", nargs='+', help="Input isea4t code(s), e.g., isea4t2geojson 131023133313201333311333 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(isea4t2geojson(args.isea4t))
    print(geojson_data)


def isea3h2geojson(isea3h_ids):
    if isinstance(isea3h_ids, str):
        isea3h_ids = [isea3h_ids]
    features = []
    if (platform.system() == 'Windows'):
        for isea3h_id in isea3h_ids:
            try:
                isea3h_dggs = Eaggr(Model.ISEA3H)
                isea3h_cell = DggsCell(isea3h_id)
                cell_polygon = isea3h_cell_to_polygon(isea3h_dggs, isea3h_cell)
                cell_centroid = cell_polygon.centroid
                center_lat = round(cell_centroid.y, 7)
                center_lon = round(cell_centroid.x, 7)
                cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]), 3)
                cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
                isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)
                cell_accuracy = isea3h2point._accuracy
                avg_edge_len = cell_perimeter / 6
                cell_resolution = isea3h_accuracy_res_dict.get(cell_accuracy)
                if (cell_resolution == 0):
                    avg_edge_len = cell_perimeter / 3
                if cell_accuracy == 0.0:
                    if round(avg_edge_len, 2) == 0.06:
                        cell_resolution = 33
                    elif round(avg_edge_len, 2) == 0.03:
                        cell_resolution = 34
                    elif round(avg_edge_len, 2) == 0.02:
                        cell_resolution = 35
                    elif round(avg_edge_len, 2) == 0.01:
                        cell_resolution = 36
                    elif round(avg_edge_len, 3) == 0.007:
                        cell_resolution = 37
                    elif round(avg_edge_len, 3) == 0.004:
                        cell_resolution = 38
                    elif round(avg_edge_len, 3) == 0.002:
                        cell_resolution = 39
                    elif round(avg_edge_len, 3) <= 0.001:
                        cell_resolution = 40
                feature = {
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                        "isea3h": isea3h_id,
                        "resolution": cell_resolution,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "avg_edge_len": round(avg_edge_len, 3),
                        "cell_area": cell_area
                    }
                }
                features.append(feature)
            except Exception:
                continue
        feature_collection = {
            "type": "FeatureCollection",
            "features": features
        }
        return feature_collection

def isea3h2geojson_cli():
    """
    Command-line interface for isea3h2geojson supporting multiple ISEA3H cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert ISEA3H ID(s) to GeoJSON")
    parser.add_argument("isea3h", nargs='+', help="Input ISEA3H cell ID(s), e.g., isea3h2geojson 1327916769,-55086 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(isea3h2geojson(args.isea3h))
    print(geojson_data)


def dggrid2geojson(dggrid_id,dggs_type,resolution):
    if (platform.system() == 'Linux'):
        dggrid_instance = DGGRIDv7(executable='/usr/local/bin/dggrid', working_dir='.', capture_logs=False, silent=True, tmp_geo_out_legacy=False, debug=False)
        dggrid_cell =  dggrid_instance.grid_cell_polygons_from_cellids([dggrid_id],dggs_type,resolution,split_dateline=True)    
      
        gdf = dggrid_cell.set_geometry("geometry")  # Ensure the geometry column is set
        # Check and set CRS to EPSG:4326 if needed
        if gdf.crs is None:
            gdf.set_crs(epsg=4326, inplace=True)
        elif not gdf.crs.equals("EPSG:4326"):
            gdf = gdf.to_crs(epsg=4326)
        
        feature_collection = gdf.to_json()
        return feature_collection


def dggrid2geojson_cli():
    """
    Command-line interface for dggrid2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert DGGRID code to GeoJSON. \
                                     Usage: dggrid2geojson <SEQNUM> <dggs_type> <res>. \
                                     Ex: dggrid2geojson 783229476878 ISEA7H 13")
    parser.add_argument("dggrid", help="Input DGGRID code in SEQNUM format")
    parser.add_argument("type", choices=dggs_types, help="Select a DGGS type from the available options.")
    parser.add_argument("res", type=int, help="resolution")
    # parser.add_argument("address", choices=input_address_types, help="Address type")

    args = parser.parse_args()
    geojson_data = dggrid2geojson(args.dggrid,args.type, args.res)
    print(geojson_data)


def ease2geojson(ease_ids):
    if isinstance(ease_ids, str):
        ease_ids = [ease_ids]
    ease_features = []
    for ease_id in ease_ids:
        try:
            level = int(ease_id[1])
            level_spec = levels_specs[level]
            n_row = level_spec["n_row"]
            n_col = level_spec["n_col"]
            geo = grid_ids_to_geos([ease_id])
            center_lon, center_lat = geo['result']['data'][0]
            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))
            cell_polygon = Polygon([
                [cell_min_lon, cell_min_lat],
                [cell_max_lon, cell_min_lat],
                [cell_max_lon, cell_max_lat],
                [cell_min_lon, cell_max_lat],
                [cell_min_lon, cell_min_lat]
            ])
            resolution = level
            num_edges = 4
            ease_feature = geodesic_dggs_to_feature("ease", ease_id, resolution, cell_polygon, num_edges)
            ease_features.append(ease_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": ease_features
    }

def ease2geojson_cli():
    """
    Command-line interface for ease2geojson supporting multiple EASE-DGGS codes.
    """
    parser = argparse.ArgumentParser(description="Convert EASE-DGGS code(s) to GeoJSON")
    parser.add_argument("ease", nargs='+', help="Input EASE-DGGS code(s), e.g., ease2geojson L4.165767.02.02.20.71 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(ease2geojson(args.ease))
    print(geojson_data)

def qtm2geojson(qtm_ids):
    """
    Convert a list of QTM cell IDs to a GeoJSON FeatureCollection.
    Accepts a single qtm_id (string) or a list of qtm_ids.
    Skips invalid or error-prone cells.
    """
    if isinstance(qtm_ids, str):
        qtm_ids = [qtm_ids]
    qtm_features = []
    for qtm_id in qtm_ids:
        try:
            facet = qtm_id_to_facet(qtm_id)
            cell_polygon = constructGeometry(facet)
            resolution = len(qtm_id)
            num_edges = 3
            qtm_feature = geodesic_dggs_to_feature("qtm", qtm_id, resolution, cell_polygon, num_edges)
            qtm_features.append(qtm_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": qtm_features
    }

def qtm2geojson_cli():
    """
    Command-line interface for qtm2geojson supporting multiple QTM cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert QTM cell ID(s) to GeoJSON")
    parser.add_argument("qtm", nargs='+', help="Input QTM cell ID(s), e.g., qtm2geojson 42012321 42012322")
    args = parser.parse_args()
    geojson_data = json.dumps(qtm2geojson(args.qtm))
    print(geojson_data)

    
def olc2geojson(olc_id):
    # Decode the Open Location Code into a CodeArea object
    coord = olc.decode(olc_id)    
    if coord:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi        
        resolution = coord.codeLength 

        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        olc_features = []
        olc_feature = graticule_dggs_to_feature("olc",olc_id,resolution,cell_polygon)   
        olc_features.append(olc_feature)

    return {
        "type": "FeatureCollection",
        "features": olc_features
    }

def olc2geojson_cli():
    """
    Command-line interface for olc2geojson supporting multiple OLC codes.
    """
    parser = argparse.ArgumentParser(description="Convert OLC/ Google Plus Codes to GeoJSON")
    parser.add_argument("olc", nargs='+', help="Input OLC(s), e.g., olc2geojson 7P28QPG4+4P7 7P28QPG4+4P8")
    args = parser.parse_args()
    geojson_data = json.dumps(olc2geojson(args.olc))
    print(geojson_data)


def geohash2geojson(geohash_ids):
    if isinstance(geohash_ids, str):
        geohash_ids = [geohash_ids]
    geohash_features = []
    for geohash_id in geohash_ids:
        try:
            bbox = geohash.bbox(geohash_id)
            if bbox:
                min_lat, min_lon = bbox['s'], bbox['w']
                max_lat, max_lon = bbox['n'], bbox['e']
                resolution = len(geohash_id)
                cell_polygon = Polygon([
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat]
                ])
                geohash_feature = graticule_dggs_to_feature("geohash", geohash_id, resolution, cell_polygon)
                geohash_features.append(geohash_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": geohash_features
    }

def geohash2geojson_cli():
    """
    Command-line interface for geohash2geojson supporting multiple Geohash cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert Geohash cell ID(s) to GeoJSON")
    parser.add_argument("geohash", nargs='+', help="Input Geohash cell ID(s), e.g., geohash2geojson w3gvk1td8 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(geohash2geojson(args.geohash))
    print(geojson_data)


def mgrs2geojson(mgrs_ids):
    if isinstance(mgrs_ids, str):
        mgrs_ids = [mgrs_ids]
    mgrs_features = []
    for mgrs_id in mgrs_ids:
        try:
            min_lat, min_lon, max_lat, max_lon, resolution = mgrs.mgrscell(mgrs_id)
            cell_polygon = Polygon([
                (min_lon, min_lat),
                (max_lon, min_lat),
                (max_lon, max_lat),
                (min_lon, max_lat),
                (min_lon, min_lat)
            ])
            mgrs_feature = graticule_dggs_to_feature("mgrs", mgrs_id, resolution, cell_polygon)
            try:
                gzd_json_path = os.path.join(os.path.dirname(__file__), '../generator/gzd.geojson')
                with open(gzd_json_path, 'r') as f:
                    gzd_data = json.load(f)
                gzd_features = gzd_data["features"]
                gzd_feature = [feature for feature in gzd_features if feature["properties"].get("gzd") == mgrs_id[:3]][0]
                gzd_geom = shape(gzd_feature["geometry"])
                if mgrs_id[2] not in {"A", "B", "Y", "Z"}:
                    if cell_polygon.intersects(gzd_geom) and not gzd_geom.contains(cell_polygon):
                        intersected_polygon = cell_polygon.intersection(gzd_geom)
                        if intersected_polygon:
                            mgrs_feature = graticule_dggs_to_feature("mgrs", mgrs_id, resolution, intersected_polygon)
            except:
                pass
            mgrs_features.append(mgrs_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": mgrs_features
    }

def mgrs2geojson_cli():
    """
    Command-line interface for mgrs2geojson supporting multiple MGRS cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert MGRS cell ID(s) to GeoJSON")
    parser.add_argument("mgrs", nargs='+', help="Input MGRS cell ID(s), e.g., mgrs2geojson 48PXS866916 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(mgrs2geojson(args.mgrs))
    print(geojson_data)


def georef2geojson(georef_ids):
    if isinstance(georef_ids, str):
        georef_ids = [georef_ids]
    georef_features = []
    for georef_id in georef_ids:
        try:
            center_lat, center_lon, min_lat, min_lon, max_lat, max_lon, resolution = georef.georefcell(georef_id)
            if center_lat:
                cell_polygon = Polygon([
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat]
                ])
                georef_feature = graticule_dggs_to_feature("georef", georef_id, resolution, cell_polygon)
                georef_features.append(georef_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": georef_features
    }

def georef2geojson_cli():
    """
    Command-line interface for georef2geojson supporting multiple GEOREF codes.
    """
    parser = argparse.ArgumentParser(description="Convert GEOREF code(s) to GeoJSON")
    parser.add_argument("georef", nargs='+', help="Input GEOREF code(s), e.g., georef2geojson VGBL42404651 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(georef2geojson(args.georef))
    print(geojson_data)


def tilecode2geojson(tilecode_ids):
    if isinstance(tilecode_ids, str):
        tilecode_ids = [tilecode_ids]
    tilecode_features = []
    for tilecode_id in tilecode_ids:
        try:
            match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)
            if not match:
                continue
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                cell_polygon = Polygon([
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat]
                ])
                resolution = z
                tilecode_feature = graticule_dggs_to_feature("tilecode_id", tilecode_id, resolution, cell_polygon)
                tilecode_features.append(tilecode_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }

def tilecode2geojson_cli():
    """
    Command-line interface for tilecode2geojson supporting multiple Tilecodes.
    """
    parser = argparse.ArgumentParser(description="Convert Tilecode(s) to GeoJSON")
    parser.add_argument("tilecode_id", nargs='+', help="Input Tilecode(s), e.g. z0x0y0 z1x1y1")
    args = parser.parse_args()
    geojson_data = json.dumps(tilecode2geojson(args.tilecode_id))
    print(geojson_data)


def quadkey2geojson(quadkey_ids):
    if isinstance(quadkey_ids, str):
        quadkey_ids = [quadkey_ids]
    quadkey_features = []
    for quadkey_id in quadkey_ids:
        try:
            tile = mercantile.quadkey_to_tile(quadkey_id)
            z = tile.z
            x = tile.x
            y = tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                cell_polygon = Polygon([
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat]
                ])
                resolution = z
                quadkey_feature = graticule_dggs_to_feature("quadkey", quadkey_id, resolution, cell_polygon)
                quadkey_features.append(quadkey_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": quadkey_features
    }

def quadkey2geojson_cli():
    """
    Command-line interface for quadkey2geojson supporting multiple Quadkeys.
    """
    parser = argparse.ArgumentParser(description="Convert Quadkey(s) to GeoJSON")
    parser.add_argument("quadkey", nargs='+', help="Input Quadkey(s), e.g. 13223011131020220011133 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(quadkey2geojson(args.quadkey))
    print(geojson_data)


def maidenhead2geojson(maidenhead_ids):
    if isinstance(maidenhead_ids, str):
        maidenhead_ids = [maidenhead_ids]
    maidenhead_features = []
    for maidenhead_id in maidenhead_ids:
        try:
            _, _, min_lat, min_lon, max_lat, max_lon, _ = maidenhead.maidenGrid(maidenhead_id)
            if min_lat:
                resolution = int(len(maidenhead_id) / 2)
                cell_polygon = Polygon([
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat]
                ])
                maidenhead_feature = graticule_dggs_to_feature("maidenhead", maidenhead_id, resolution, cell_polygon)
                maidenhead_features.append(maidenhead_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": maidenhead_features
    }

def maidenhead2geojson_cli():
    """
    Command-line interface for maidenhead2geojson supporting multiple Maidenhead cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert Maidenhead cell ID(s) to GeoJSON")
    parser.add_argument("maidenhead", nargs='+', help="Input Maidenhead cell ID(s), e.g., maidenhead2geojson OK30is46 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(maidenhead2geojson(args.maidenhead))
    print(geojson_data)


def gars2geojson(gars_ids):
    if isinstance(gars_ids, str):
        gars_ids = [gars_ids]
    gars_features = []
    for gars_id in gars_ids:
        try:
            gars_grid = garsgrid.GARSGrid(gars_id)
            wkt_polygon = gars_grid.polygon
            if wkt_polygon:
                resolution_minute = gars_grid.resolution
                resolution = 1
                if resolution_minute == 30:
                    resolution = 1
                elif resolution_minute == 15:
                    resolution = 2
                elif resolution_minute == 5:
                    resolution = 3
                elif resolution_minute == 1:
                    resolution = 4
                cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
                gars_feature = graticule_dggs_to_feature("gars", gars_id, resolution, cell_polygon)
                gars_features.append(gars_feature)
        except Exception:
            continue
    return {
        "type": "FeatureCollection",
        "features": gars_features
    }

def gars2geojson_cli():
    """
    Command-line interface for gars2geojson supporting multiple GARS cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert GARS cell ID(s) to GeoJSON")
    parser.add_argument("gars", nargs='+', help="Input GARS cell ID(s), e.g., gars2geojson 574JK1918 ...")
    args = parser.parse_args()
    geojson_data = json.dumps(gars2geojson(args.gars))
    print(geojson_data)
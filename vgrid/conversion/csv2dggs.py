import os,argparse,json, re
import pandas as pd
from tqdm import tqdm
import h3
from urllib.parse import urlparse
import requests
from io import StringIO

from shapely.geometry import Polygon, mapping
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import  chunk_size, geodesic_dggs_to_feature

from vgrid.utils import s2, olc, geohash, georef, mgrs, mercantile, maidenhead
from vgrid.utils.gars import garsgrid

from shapely.wkt import loads

from vgrid.utils.antimeridian import fix_polygon

from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.conversion.dggs2geojson import rhealpix_cell_to_polygon

import platform
if (platform.system() == 'Windows'):   
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.generator.isea4tgrid import fix_isea4t_wkt, fix_isea4t_antimeridian_cells
    
    from vgrid.conversion.dggs2geojson import isea3h_cell_to_polygon

from vgrid.utils.easedggs.constants import levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos

from vgrid.generator.settings import isea3h_accuracy_res_dict, geodesic_dggs_to_feature, graticule_dggs_to_feature
from vgrid.utils.qtm import constructGeometry, qtm_id_to_facet


from pyproj import Geod
geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID


#################################################################################
#  H3
#################################################################################
def h32feature(h3_id):
    """Convert H3 cell ID to a GeoJSON Polygon."""
    cell_boundary = h3.cell_to_boundary(h3_id)
    if cell_boundary:
        filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
        # Reverse lat/lon to lon/lat for GeoJSON compatibility
        reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
        cell_polygon = Polygon(reversed_boundary)
        resolution = h3.get_resolution(h3_id)             
        num_edges = 6
        if h3.is_pentagon(h3_id):
            num_edges = 6         
        h3_feature = geodesic_dggs_to_feature("h3",h3_id,resolution,cell_polygon,num_edges)              
        
        return h3_feature
    
def csv2h3(csv_file, column='h3'):
    """Convert CSV file containing H3 cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            h3_id = row[column]
            h3_feature = h32feature(h3_id)
            if h3_feature:
                geojson_features.append(h3_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid H3 cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2h3_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with H3 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with H3 column")
    parser.add_argument("id", help="Name of the H3 column (default: 'h3')", default='h3')
    args = parser.parse_args()
    h3_csv = args.csv
    h3_id = args.id

    h3_geojson = csv2h3(h3_csv, h3_id)
    geojson_name = os.path.splitext(os.path.basename(h3_csv))[0]
    geojson_path = f"{geojson_name}2h3.geojson"

    with open(geojson_path, "w") as f:
        json.dump(h3_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")

#################################################################################
#  S2
#################################################################################
def s22feature(s2_token):
    # Create an S2 cell from the given cell ID
    cell_id = s2.CellId.from_token(s2_token)
    cell = s2.Cell(cell_id)
    if cell:
        # Get the vertices of the cell (4 vertices for a rectangular cell)
        vertices = [cell.get_vertex(i) for i in range(4)]
        # Prepare vertices in (longitude, latitude) format for Shapely
        shapely_vertices = []
        for vertex in vertices:
            lat_lng = s2.LatLng.from_point(vertex)  # Convert Point to LatLng
            longitude = lat_lng.lng().degrees  # Access longitude in degrees
            latitude = lat_lng.lat().degrees   # Access latitude in degrees
            shapely_vertices.append((longitude, latitude))

        # Close the polygon by adding the first vertex again
        shapely_vertices.append(shapely_vertices[0])  # Closing the polygon
        # Create a Shapely Polygon
        cell_polygon = fix_polygon(Polygon(shapely_vertices)) # Fix antimeridian
        resolution = cell_id.level()            
        num_edges = 4        
        
        s2_feature = geodesic_dggs_to_feature("s2",s2_token,resolution,cell_polygon,num_edges)              
        return s2_feature

def is_url(path):
    """Check if the given path is a URL."""
    try:
        result = urlparse(path)
        return all([result.scheme, result.netloc])
    except:
        return False

def read_csv_file(csv_file):
    """Read CSV file from either local path or URL."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()  # Raise an exception for bad status codes
            return pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        return pd.read_csv(csv_file)

def csv2s2(csv_file, column='s2'):
    """Convert CSV file containing S2 cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            s2_id = row[column]
            s2_feature = s22feature(s2_id)
            if s2_feature:
                geojson_features.append(s2_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid S2 cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2s2_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with S2 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with S2 column")
    parser.add_argument("id", help="Name of the S2 column (default: 's2')", default='s2')
    args = parser.parse_args()
    s2_csv = args.csv
    s2_id = args.id

    s2_geojson = csv2s2(s2_csv, s2_id)
    geojson_name = os.path.splitext(os.path.basename(s2_csv))[0]
    geojson_path = f"{geojson_name}2s2.geojson"

    with open(geojson_path, "w") as f:
        json.dump(s2_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")

#################################################################################
#  Rhealpix
#################################################################################
def rhealpix2feature(rhealpix_id):
    rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
    rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3) 
    rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
    if rhealpix_cell:
        resolution = rhealpix_cell.resolution        
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
        num_edges = 4
        if rhealpix_cell.ellipsoidal_shape() == 'dart':
            num_edges = 3
        rhealpix_feature = geodesic_dggs_to_feature("rhealpix",rhealpix_id,resolution,cell_polygon,num_edges)                
        return rhealpix_feature

def csv2rhealpix(csv_file, column='rhealpix'):
    """Convert CSV file containing rHEALPix cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            rhealpix_id = row[column]
            rhealpix_feature = rhealpix2feature(rhealpix_id)
            if rhealpix_feature:
                geojson_features.append(rhealpix_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid rHEALPix cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2rhealpix_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with rhealpix column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with rhealpix column")
    parser.add_argument("id", help="Name of the rhealpix column (default: 'rhealpix')", default='rhealpix')
    args = parser.parse_args()
    rhealpix_csv = args.csv
    rhealpix_id = args.id

    rhealpix_geojson = csv2rhealpix(rhealpix_csv, rhealpix_id)
    geojson_name = os.path.splitext(os.path.basename(rhealpix_csv))[0]
    geojson_path = f"{geojson_name}2rhealpix.geojson"

    with open(geojson_path, "w") as f:
        json.dump(rhealpix_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")

#################################################################################
#  Open-Eaggr ISEA4T
#################################################################################
def isea4t2feature(isea4t_id):
    if (platform.system() == 'Windows'): 
        isea4t_dggs = Eaggr(Model.ISEA4T)
        cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(isea4t_id),ShapeStringFormat.WKT)
        cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
        if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14')\
            or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
            cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
        
        if cell_to_shape_fixed:
            resolution = len(isea4t_id)-2
            num_edges = 3
            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
            isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,resolution,cell_polygon,num_edges)   
            return isea4t_feature

def csv2isea4t(csv_file, column='isea4t'):
    """Convert CSV file containing ISEA4T cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            isea4t_id = row[column]
            isea4t_feature = isea4t2feature(isea4t_id)
            if isea4t_feature:
                geojson_features.append(isea4t_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid ISEA4T cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2isea4t_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with ISEA4T column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with ISEA4T column")
    parser.add_argument("id", help="Name of the ISEA4T column (default: 'isea4t')", default='isea4t')
    args = parser.parse_args()
    isea4t_csv = args.csv
    isea4t_id = args.id

    isea4t_geojson = csv2isea4t(isea4t_csv, isea4t_id)
    geojson_name = os.path.splitext(os.path.basename(isea4t_csv))[0]
    geojson_path = f"{geojson_name}2isea4t.geojson"

    with open(geojson_path, "w") as f:
        json.dump(isea4t_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")

#################################################################################
#  Open-Eaggr ISEA3H
#################################################################################
def isea3h2feature(isea3h_id):
    if (platform.system() == 'Windows'): 
        isea3h_dggs = Eaggr(Model.ISEA3H)
        cell_polygon = isea3h_cell_to_polygon(isea3h_id)
        if cell_polygon:
    
            cell_centroid = cell_polygon.centroid
            center_lat =  round(cell_centroid.y, 7)
            center_lon = round(cell_centroid.x, 7)
            
            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            
            isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(isea3h_id))      
            accuracy = isea3h2point._accuracy
                
            avg_edge_len = cell_perimeter / 6
            resolution  = isea3h_accuracy_res_dict.get(accuracy)
            
            if (resolution == 0): # icosahedron faces at resolution = 0
                avg_edge_len = cell_perimeter / 3
            
            if accuracy == 0.0:
                if round(avg_edge_len,2) == 0.06:
                    resolution = 33
                elif round(avg_edge_len,2) == 0.03:
                    resolution = 34
                elif round(avg_edge_len,2) == 0.02:
                    resolution = 35
                elif round(avg_edge_len,2) == 0.01:
                    resolution = 36
                
                elif round(avg_edge_len,3) == 0.007:
                    resolution = 37
                elif round(avg_edge_len,3) == 0.004:
                    resolution = 38
                elif round(avg_edge_len,3) == 0.002:
                    resolution = 39
                elif round(avg_edge_len,3) <= 0.001:
                    resolution = 40
                    
            isea3h_feature = {
                "type": "Feature",
                "geometry": mapping(cell_polygon),
                "properties": {
                        "isea3h": isea3h_id,
                        "resolution": resolution,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "avg_edge_len": round(avg_edge_len,3),
                        "cell_area": cell_area
                        }
            }
            return isea3h_feature

      
def csv2isea3h(csv_file, id_col=None):
    if not os.path.exists(csv_file):
        print(f"Error: Input file {csv_file} does not exist.")
        return    
    
    if id_col is None:
        id_col = 'isea3h'
        
    try:
        df = pd.read_csv(csv_file, dtype=str)  # Read entire file
        
        if id_col not in df.columns:
            print(f"Error: Column '{id_col}' is missing in the input CSV. Please check and try again.")
            return
        
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            isea3h_id = row[id_col]
            isea3h_feature = isea3h2feature(isea3h_id)
            if isea3h_feature:
                isea3h_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                geojson_features.append(isea3h_feature)
        except Exception as e:
            print(f" Skipping row {row.to_dict()}: {e}")
    
    isea3h_geojson = {"type": "FeatureCollection", "features": geojson_features}
    return isea3h_geojson

def csv2isea3h_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with ISEA3H column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with ISEA3H column")
    parser.add_argument("id", help="Name of the ISEA3H column (default: 'isea3h')", default='isea3h')
    args = parser.parse_args()
    isea3h_csv = args.csv
    isea3h_id = args.id

    isea3h_geojson = csv2isea3h(isea3h_csv, isea3h_id)
    geojson_name = os.path.splitext(os.path.basename(isea3h_csv))[0]
    geojson_path = f"{geojson_name}2isea3h.geojson"

    with open(geojson_path, "w") as f:
        json.dump(isea3h_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")
    
#################################################################################
#  EASE-DGGS
#################################################################################
def ease2feature(ease_id):
    level = int(ease_id[1])  # Get the level (e.g., 'L0' -> 0)
    # Get level specs
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

    if cell_polygon:
        resolution = level
        num_edges = 4
        ease_feature = geodesic_dggs_to_feature("ease",ease_id,resolution,cell_polygon,num_edges)   
        return ease_feature

def csv2ease(csv_file, column='ease'):
    """Convert CSV file containing EASE cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            ease_id = row[column]
            ease_feature = ease2feature(ease_id)
            if ease_feature:
                geojson_features.append(ease_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid EASE cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2ease_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with EASE column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with EASE column")
    parser.add_argument("id", help="Name of the EASE column (default: 'ease')", default='ease')
    args = parser.parse_args()
    ease_csv = args.csv
    ease_id = args.id

    ease_geojson = csv2ease(ease_csv, ease_id)
    geojson_name = os.path.splitext(os.path.basename(ease_csv))[0]
    geojson_path = f"{geojson_name}2ease.geojson"

    with open(geojson_path, "w") as f:
        json.dump(ease_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  QTM
#################################################################################
def qtm2feature(qtm_id):
    facet = qtm_id_to_facet(qtm_id)
    cell_polygon = constructGeometry(facet)    
    if cell_polygon:
        resolution = len(qtm_id)
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id,resolution,cell_polygon,num_edges)   
        return qtm_feature
    
def csv2qtm(csv_file, column='qtm'):
    """Convert CSV file containing QTM cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            qtm_id = row[column]
            qtm_feature = qtm2feature(qtm_id)
            if qtm_feature:
                geojson_features.append(qtm_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid QTM cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2qtm_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with qtm column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with qtm column")
    parser.add_argument("id", help="Name of the qtm column (default: 'qtm')", default='qtm')
    args = parser.parse_args()
    qtm_csv = args.csv
    qtm_id = args.id

    qtm_geojson = csv2qtm(qtm_csv, qtm_id)
    geojson_name = os.path.splitext(os.path.basename(qtm_csv))[0]
    geojson_path = f"{geojson_name}2qtm.geojson"

    with open(geojson_path, "w") as f:
        json.dump(qtm_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")



#################################################################################
#  OLC
#################################################################################
def olc2feature(olc_id):
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
        olc_feature = graticule_dggs_to_feature("olc",olc_id,resolution,cell_polygon)   
        return olc_feature
    
def csv2olc(csv_file, column='olc'):
    """Convert CSV file containing OLC cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            olc_id = row[column]
            olc_feature = olc2feature(olc_id)
            if olc_feature:
                geojson_features.append(olc_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid OLC cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2olc_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with OLC column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with OLC column")
    parser.add_argument("id", help="Name of the OLC column (default: 'olc')", default='olc')
    args = parser.parse_args()
    olc_csv = args.csv
    olc_id = args.id

    olc_geojson = csv2olc(olc_csv, olc_id)
    geojson_name = os.path.splitext(os.path.basename(olc_csv))[0]
    geojson_path = f"{geojson_name}2olc.geojson"

    with open(geojson_path, "w") as f:
        json.dump(olc_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  Geohash
#################################################################################
def geohash2feature(geohash_id):
    bbox =  geohash.bbox(geohash_id)
    if bbox:
        min_lat, min_lon = bbox['s'], bbox['w']  # Southwest corner
        max_lat, max_lon = bbox['n'], bbox['e']  # Northeast corner        
        resolution =  len(geohash_id)

        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        geohash_feature = graticule_dggs_to_feature("geohash",geohash_id,resolution,cell_polygon)   
        return geohash_feature
    
def csv2geohash(csv_file, column='geohash'):
    """Convert CSV file containing Geohash cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            geohash_id = row[column]
            geohash_feature = geohash2feature(geohash_id)
            if geohash_feature:
                geojson_features.append(geohash_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid Geohash cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2geohash_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with Geohash column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with Geohash column")
    parser.add_argument("id", help="Name of the Geohash column (default: 'geohash')", default='geohash')
    args = parser.parse_args()
    geohash_csv = args.csv
    geohash_id = args.id

    geohash_geojson = csv2geohash(geohash_csv, geohash_id)
    geojson_name = os.path.splitext(os.path.basename(geohash_csv))[0]
    geojson_path = f"{geojson_name}2geohash.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geohash_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  GEOREF
#################################################################################
def georef2feature(georef_id):
    # Need to check georef.georefcell(georef_id) function
    center_lat, center_lon, min_lat, min_lon, max_lat, max_lon,resolution = georef.georefcell(georef_id)
    if center_lat:
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        georef_feature = graticule_dggs_to_feature("georef",georef_id,resolution,cell_polygon) 
        return   georef_feature
    
def csv2georef(csv_file, column='georef'):
    """Convert CSV file containing Georef cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            georef_id = row[column]
            georef_feature = georef2feature(georef_id)
            if georef_feature:
                geojson_features.append(georef_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid Georef cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2georef_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with GEOREF column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with GEOREF column")
    parser.add_argument("id", help="Name of the GEOREF column (default: 'georef')", default='georef')
    args = parser.parse_args()
    georef_csv = args.csv
    georef_id = args.id

    georef_geojson = csv2georef(georef_csv, georef_id)
    geojson_name = os.path.splitext(os.path.basename(georef_csv))[0]
    geojson_path = f"{geojson_name}2georef.geojson"

    with open(geojson_path, "w") as f:
        json.dump(georef_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")
    
    
#################################################################################
#  MGRS
#################################################################################
def mgrs2feature(mgrs_id):
    # Need to check if MGRS cell is intersectd by GZD
    min_lat, min_lon, max_lat, max_lon,resolution = mgrs.mgrscell(mgrs_id)
    if min_lat:
        # Define the polygon based on the bounding box       
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])

        mgrs_feature = graticule_dggs_to_feature("georef",mgrs_id,resolution,cell_polygon) 
        return   mgrs_feature
        
    
def csv2mgrs(csv_file, column='mgrs'):
    """Convert CSV file containing MGRS cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            mgrs_id = row[column]
            mgrs_feature = mgrs2feature(mgrs_id)
            if mgrs_feature:
                geojson_features.append(mgrs_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid MGRS cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2mgrs_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with mgrs column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with mgrs column")
    parser.add_argument("id", help="Name of the mgrs column (default: 'mgrs')", default='mgrs')
    args = parser.parse_args()
    mgrs_csv = args.csv
    mgrs_id = args.id

    mgrs_geojson = csv2mgrs(mgrs_csv, mgrs_id)
    geojson_name = os.path.splitext(os.path.basename(mgrs_csv))[0]
    geojson_path = f"{geojson_name}2mgrs.geojson"

    with open(geojson_path, "w") as f:
        json.dump(mgrs_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")

#################################################################################
#  Tilecode
#################################################################################
def tilecode2feature(tilecode_id):
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)
    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile in (west, south, east, north)
    bounds = mercantile.bounds(x, y, z)    
    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        if cell_polygon:
            resolution = z
            tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id,resolution,cell_polygon)           
            return tilecode_feature
        
def csv2tilecode(csv_file, column='tilecode'):
    """Convert CSV file containing Tilecode cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            tilecode_id = row[column]
            tilecode_feature = tilecode2feature(tilecode_id)
            if tilecode_feature:
                geojson_features.append(tilecode_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid Tilecode cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2tilecode_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with tilecode column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with tilecode column")
    parser.add_argument("id", help="Name of the tilecode column (default: 'tilecode')", default='tilecode')
    args = parser.parse_args()
    tilecode_csv = args.csv
    tilecode_id = args.id

    tilecode_geojson = csv2tilecode(tilecode_csv, tilecode_id)
    geojson_name = os.path.splitext(os.path.basename(tilecode_csv))[0]
    geojson_path = f"{geojson_name}2tilecode.geojson"

    with open(geojson_path, "w") as f:
        json.dump(tilecode_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  Quadkey
#################################################################################
def quadkey2feature(quadkey_id):
    tile = mercantile.quadkey_to_tile(quadkey_id)    
    # Format as tilecode
    z = tile.z
    x = tile.x
    y = tile.y
    # Get the bounds of the tile in (west, south, east, north)
    bounds = mercantile.bounds(x, y, z)    
    
    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east

        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        resolution = z
        quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id,resolution,cell_polygon)   
        return quadkey_feature
        
def csv2quadkey(csv_file, column='quadkey'):
    """Convert CSV file containing Quadkey cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            quadkey_id = row[column]
            quadkey_feature = quadkey2feature(quadkey_id)
            if quadkey_feature:
                geojson_features.append(quadkey_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid Quadkey cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2quadkey_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with quadkey column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with quadkey column")
    parser.add_argument("id", help="Name of the quadkey column (default: 'quadkey')", default='quadkey')
    args = parser.parse_args()
    quadkey_csv = args.csv
    quadkey_id = args.id

    quadkey_geojson = csv2quadkey(quadkey_csv, quadkey_id)
    geojson_name = os.path.splitext(os.path.basename(quadkey_csv))[0]
    geojson_path = f"{geojson_name}2quadkey.geojson"

    with open(geojson_path, "w") as f:
        json.dump(quadkey_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  Maidenhead
#################################################################################
def maidenhead2feature(maidenhead_id):
    # Decode the Open Location Code into a CodeArea object
    _, _, min_lat, min_lon, max_lat, max_lon, _ = maidenhead.maidenGrid(maidenhead_id)
    if min_lat:
        resolution = int(len(maidenhead_id)/2)   
        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        maidenhead_feature = graticule_dggs_to_feature("gars",maidenhead_id,resolution,cell_polygon) 
        return maidenhead_feature
    
def csv2maidenhead(csv_file, column='maidenhead'):
    """Convert CSV file containing Maidenhead cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            maidenhead_id = row[column]
            maidenhead_feature = maidenhead2feature(maidenhead_id)
            if maidenhead_feature:
                geojson_features.append(maidenhead_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid Maidenhead cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2maidenhead_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with Maidenhead column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with Maidenhead column")
    parser.add_argument("id", help="Name of the Maidenhead column (default: 'maidenhead')", default='maidenhead')
    args = parser.parse_args()
    maidenhead_csv = args.csv
    maidenhead_id = args.id

    maidenhead_geojson = csv2maidenhead(maidenhead_csv, maidenhead_id)
    geojson_name = os.path.splitext(os.path.basename(maidenhead_csv))[0]
    geojson_path = f"{geojson_name}2maidenhead.geojson"

    with open(geojson_path, "w") as f:
        json.dump(maidenhead_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  GARS
#################################################################################
def gars2feature(gars_id):
    gars_grid = garsgrid.GARSGrid(gars_id)
    wkt_polygon = gars_grid.polygon
    
    if wkt_polygon:
        # Convert minute-based resolution to 1-4 scale
        resolution_minute = gars_grid.resolution
        if resolution_minute == 30:
            resolution = 1
        elif resolution_minute == 15:
            resolution = 2
        elif resolution_minute == 5:
            resolution = 3
        elif resolution_minute == 1:
            resolution = 4
        else:
            resolution = 1  # Default to level 1 if unknown
            
        cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
        gars_feature = graticule_dggs_to_feature("gars",gars_id,resolution,cell_polygon)   
        return gars_feature
    
    
def csv2gars(csv_file, column='gars'):
    """Convert CSV file containing GARS cell IDs to GeoJSON."""
    if is_url(csv_file):
        try:
            response = requests.get(csv_file)
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text))
        except requests.RequestException as e:
            print(f"Error: Failed to download CSV from URL {csv_file}: {str(e)}")
            return None
    else:
        if not os.path.exists(csv_file):
            print(f"Error: Input file {csv_file} does not exist.")
            return None
        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return None
    
    if column not in df.columns:
        print(f"Error: Column '{column}' is missing in the input CSV. Please check and try again.")
        return None
    
    geojson_features = []    
    for _, row in df.iterrows():
        try:
            gars_id = row[column]
            gars_feature = gars2feature(gars_id)
            if gars_feature:
                geojson_features.append(gars_feature)
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue
    
    if not geojson_features:
        print("No valid GARS cells found in the CSV file.")
        return None
    
    geojson = {
        'type': 'FeatureCollection',
        'features': geojson_features
    }
    return geojson

def csv2gars_cli():
    parser = argparse.ArgumentParser(description="Convert CSV with GARS column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with GARS column")
    parser.add_argument("id", help="Name of the GARS column (default: 'gars')", default='gars')
    args = parser.parse_args()
    gars_csv = args.csv
    gars_id = args.id

    gars_geojson = csv2gars(gars_csv, gars_id)
    geojson_name = os.path.splitext(os.path.basename(gars_csv))[0]
    geojson_path = f"{geojson_name}2gars.geojson"

    with open(geojson_path, "w") as f:
        json.dump(gars_geojson, f)
    
    print(f"GeoJSON saved to {geojson_path}")
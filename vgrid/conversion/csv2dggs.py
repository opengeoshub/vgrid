import os,argparse,json, re
import pandas as pd
from tqdm import tqdm
import h3

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
    
def csv2h3():
    parser = argparse.ArgumentParser(description="Convert CSV with H3 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'h3' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "h3" not in columns_lower:
            print("Error: Column 'h3' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        h3_col = columns_lower["h3"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"h3": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                h3_id = row[h3_col]
                h3_feature = h32feature(h3_id)
                if h3_feature:
                    h3_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(h3_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2h3.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
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

def csv2s2():
    parser = argparse.ArgumentParser(description="Convert CSV with S2 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 's2' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "s2" not in columns_lower:
            print("Error: Column 's2' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        s2_col = columns_lower["s2"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"s2": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                s2_id = row[s2_col]
                s2_feature = s22feature(s2_id)
                if s2_feature:
                    s2_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(s2_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2s2.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
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

def csv2rhealpix():
    parser = argparse.ArgumentParser(description="Convert CSV with rhealpix column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'rhealpix' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "rhealpix" not in columns_lower:
            print("Error: Column 'rhealpix' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        rhealpix_col = columns_lower["rhealpix"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"rhealpix": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                rhealpix_id = row[rhealpix_col]
                rhealpix_feature = rhealpix2feature(rhealpix_id)
                if rhealpix_feature:
                    rhealpix_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(rhealpix_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2rhealpix.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
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

def csv2isea4t():
    parser = argparse.ArgumentParser(description="Convert CSV with ISEA4T column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'isea4t' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "isea4t" not in columns_lower:
            print("Error: Column 'isea4t' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        rhealpix_col = columns_lower["isea4t"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"isea4t": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                isea4t_id = row[rhealpix_col]
                isea4t_feature = isea4t2feature(isea4t_id)
                if isea4t_feature:
                    isea4t_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(isea4t_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2isea4t.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
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

      
def csv2isea3h():
    parser = argparse.ArgumentParser(description="Convert CSV with ISEA3H column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'isea3h' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "isea3h" not in columns_lower:
            print("Error: Column 'isea3h' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        isea3h_col = columns_lower["isea3h"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                isea3h_id = row[isea3h_col]  # Use the correct column name
                isea3h_feature = isea3h2feature(isea3h_id)
                if isea3h_feature:
                    isea3h_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(isea3h_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2isea3h.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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

def csv2ease():
    parser = argparse.ArgumentParser(description="Convert CSV with EASE column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'ease' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "ease" not in columns_lower:
            print("Error: Column 'ease' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        ease_col = columns_lower["ease"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                ease_id = row[ease_col]  # Use the correct column name
                ease_feature = ease2feature(ease_id)
                if ease_feature:
                    ease_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(ease_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2ease.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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
    
def csv2qtm():
    parser = argparse.ArgumentParser(description="Convert CSV with qtm column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'qtm' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "qtm" not in columns_lower:
            print("Error: Column 'qtm' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        rhealpix_col = columns_lower["qtm"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"qtm": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                qtm_id = row[rhealpix_col]
                qtm_feature = qtm2feature(qtm_id)
                if qtm_feature:
                    qtm_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(qtm_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2qtm.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
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
    
def csv2olc():
    parser = argparse.ArgumentParser(description="Convert CSV with OLC column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'olc' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "olc" not in columns_lower:
            print("Error: Column 'olc' (case insensitive) is missing in the input CSV. Plolc check and try again.")
            return
        
        olc_col = columns_lower["olc"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                olc_id = row[olc_col]  # Use the correct column name
                olc_feature = olc2feature(olc_id)
                if olc_feature:
                    olc_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(olc_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2olc.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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
    
def csv2geohash():
    parser = argparse.ArgumentParser(description="Convert CSV with Geohash column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'geohash' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "geohash" not in columns_lower:
            print("Error: Column 'geohash' (case insensitive) is missing in the input CSV. Plgeohash check and try again.")
            return
        
        geohash_col = columns_lower["geohash"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                geohash_id = row[geohash_col]  # Use the correct column name
                geohash_feature = geohash2feature(geohash_id)
                if geohash_feature:
                    geohash_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(geohash_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2geohash.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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
    
def csv2georef():
    parser = argparse.ArgumentParser(description="Convert CSV with GEOREF column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'georef' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "georef" not in columns_lower:
            print("Error: Column 'georef' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        georef_col = columns_lower["georef"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                georef_id = row[georef_col]  # Use the correct column name
                georef_feature = georef2feature(georef_id)
                if georef_feature:
                    georef_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(georef_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2georef.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"GeoJSON saved to {geojson_path}")
    
    
#################################################################################
#  MGRS
#################################################################################
def mgrs2feature(mgrs_id):
    # Need to check if MGRS cell is intersectd by GZD
    origin_lat, origin_lon, min_lat, min_lon, max_lat, max_lon,resolution = mgrs.mgrscell(mgrs_id)
    if origin_lat:
        # Define the polygon based on the bounding box
        origin_lat = round(origin_lat,7)
        origin_lon = round(origin_lon,7)
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])

        mgrs_feature = graticule_dggs_to_feature("georef",mgrs_id,resolution,cell_polygon) 
        return   mgrs_feature
        
    
def csv2mgrs():
    parser = argparse.ArgumentParser(description="Convert CSV with mgrs column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'mgrs' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "mgrs" not in columns_lower:
            print("Error: Column 'mgrs' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        mgrs_col = columns_lower["mgrs"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                mgrs_id = row[mgrs_col]  # Use the correct column name
                mgrs_feature = mgrs2feature(mgrs_id)
                if mgrs_feature:
                    mgrs_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(mgrs_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2mgrs.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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
        
def csv2tilecode():
    parser = argparse.ArgumentParser(description="Convert CSV with tilecode column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'tilecode' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "tilecode" not in columns_lower:
            print("Error: Column 'tilecode' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        tilecode_col = columns_lower["tilecode"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                tilecode_id = row[tilecode_col]  # Use the correct column name
                tilecode_feature = tilecode2feature(tilecode_id)
                if tilecode_feature:
                    tilecode_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(tilecode_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2tilecode.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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
        
def csv2quadkey():
    parser = argparse.ArgumentParser(description="Convert CSV with quadkey column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'quadkey' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "quadkey" not in columns_lower:
            print("Error: Column 'quadkey' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        quadkey_col = columns_lower["quadkey"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                quadkey_id = row[quadkey_col]  # Use the correct column name
                quadkey_feature = quadkey2feature(quadkey_id)
                if quadkey_feature:
                    quadkey_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(quadkey_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2quadkey.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

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
    
def csv2maidenhead():
    parser = argparse.ArgumentParser(description="Convert CSV with Maidenhead column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'maidenhead' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "maidenhead" not in columns_lower:
            print("Error: Column 'maidenhead' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        maidenhead_col = columns_lower["maidenhead"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                maidenhead_id = row[maidenhead_col]  # Use the correct column name
                maidenhead_feature = maidenhead2feature(maidenhead_id)
                if maidenhead_feature:
                    maidenhead_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(maidenhead_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2maidenhead.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"GeoJSON saved to {geojson_path}")


#################################################################################
#  GARS
#################################################################################
def gars2feature(gars_id):
    gars_grid = garsgrid.GARSGrid(gars_id)
    wkt_polygon = gars_grid.polygon
    
    if wkt_polygon:
        # # Create the bounding box coordinates for the polygon
        resolution_minute = gars_grid.resolution        
        cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
        gars_feature = graticule_dggs_to_feature("gars",gars_id,resolution_minute,cell_polygon)   
        return gars_feature
    
    
def csv2gars():
    parser = argparse.ArgumentParser(description="Convert CSV with GARS column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'gars' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    try:
        first_chunk = pd.read_csv(csv, dtype=str, nrows=1)  # Read first row to check columns
        columns_lower = {col.lower(): col for col in first_chunk.columns}  # Create a case-insensitive mapping
        
        if "gars" not in columns_lower:
            print("Error: Column 'gars' (case insensitive) is missing in the input CSV. Please check and try again.")
            return
        
        gars_col = columns_lower["gars"]  # Get the actual column name
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    geojson_features = []
    for chunk in pd.read_csv(csv, dtype=str, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                gars_id = row[gars_col]  # Use the correct column name
                gars_feature = gars2feature(gars_id)
                if gars_feature:
                    gars_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(gars_feature)
            except Exception as e:
                print(f" Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2gars.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"GeoJSON saved to {geojson_path}")
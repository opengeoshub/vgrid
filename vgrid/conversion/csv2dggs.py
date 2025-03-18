import os, argparse,json
import pandas as pd
from tqdm import tqdm
import h3

from shapely.geometry import Polygon, mapping
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import  chunk_size, geodesic_dggs_to_feature

from vgrid.utils import s2
from shapely.wkt import loads

from vgrid.generator.settings import geodesic_dggs_to_feature
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
    from vgrid.generator.settings import isea3h_accuracy_res_dict

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
                print(f"Skipping row {row.to_dict()}: {e}")
    
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
                print(f"Skipping row {row.to_dict()}: {e}")
    
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
                print(f"Skipping row {row.to_dict()}: {e}")
    
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
                print(f"Skipping row {row.to_dict()}: {e}")
    
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
                print(f"Skipping row {row.to_dict()}: {e}")

    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2isea3h.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)

    print(f"GeoJSON saved to {geojson_path}")
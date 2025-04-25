from vgrid.utils import s2, olc, geohash,  mercantile, tilecode
from vgrid.utils import qtm
import h3

from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
import platform

if (platform.system() == 'Windows'):   
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.generator.isea4tgrid import fix_isea4t_wkt, fix_isea4t_antimeridian_cells


if (platform.system() == 'Linux'):
    from vgrid.utils.dggrid4py import DGGRIDv7, dggs_types
    from vgrid.utils.dggrid4py.dggrid_runner import input_address_types


from vgrid.utils.easedggs.constants import levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos

from shapely.wkt import loads
from shapely.geometry import shape, Polygon,mapping

import json, re,os,argparse
from vgrid.generator.h3grid import fix_h3_antimeridian_cells

from vgrid.utils.antimeridian import fix_polygon

from vgrid.generator.settings import graticule_dggs_to_feature, geodesic_dggs_to_feature,isea3h_accuracy_res_dict
from vgrid.conversion.dggs2geojson import rhealpix_cell_to_polygon
from vgrid.utils.easedggs.dggs.hierarchy import _parent_to_children
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos
from vgrid.generator.geohashgrid import geohash_to_polygon

from pyproj import Geod
geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID
from collections import defaultdict

from tqdm import tqdm

#################
# H3
#################
def h3compact(geojson_data):
    # h3_ids = list(set(feature["properties"]["h3"] for feature in geojson_data.get("features", []) if "h3" in feature.get("properties", {})))
    h3_ids = [feature["properties"]["h3"] for feature in geojson_data.get("features", []) if "h3" in feature.get("properties", {})]
    h3_ids_compact = h3.compact_cells(h3_ids)
    h3_features = [] 
    for h3_id_compact in tqdm(h3_ids_compact, desc="Compacting cells "):  
        cell_boundary = h3.cell_to_boundary(h3_id_compact)   
        if cell_boundary:
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            resolution = h3.get_resolution(h3_id_compact)
            num_edges = 6
            if (h3.is_pentagon(h3_id_compact)):
                num_edges = 5
            h3_feature = geodesic_dggs_to_feature("h3",h3_id_compact,resolution,cell_polygon,num_edges)   
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features
    }
    
def h3compact_cli():
    """
    Command-line interface for h3compact.
    """
    parser = argparse.ArgumentParser(description="Compact H3 in a GeoJSON file containing a property named 'h3'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input H3 in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = h3compact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "h3_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")

def h3expand(geojson_data,resolution):
    h3_ids = [feature["properties"]["h3"] for feature in geojson_data.get("features", []) if "h3" in feature.get("properties", {})]
    h3_ids_expand = h3.uncompact_cells(h3_ids, resolution)
    h3_features = [] 
    for h3_id_expand in tqdm(h3_ids_expand, desc="Expanding cells "):
        cell_boundary = h3.cell_to_boundary(h3_id_expand)   
        if cell_boundary:
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            num_edges = 6
            if (h3.is_pentagon(h3_id_expand)):
                num_edges = 5
            h3_feature = geodesic_dggs_to_feature("h3",h3_id_expand,resolution,cell_polygon,num_edges)   
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features
    }
    
def h3expand_cli():
    """
    Command-line interface for h3expand.
    """
    parser = argparse.ArgumentParser(description="expand H3 in a GeoJSON file containing a property named 'h3'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of H3 to be expanded [0..15]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input H3 in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution < 0 or resolution > 15:
        print(f"Please select a resolution in [0..15] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = h3expand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"h3_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")

#################
# S2
#################
def s2compact(geojson_data):
    s2_tokens = [feature["properties"]["s2"] for feature in geojson_data.get("features", []) if "s2" in feature.get("properties", {})]
    s2_ids = [s2.CellId.from_token(token) for token in s2_tokens]
    if s2_ids:
        covering = s2.CellUnion(s2_ids)
        covering.normalize()
        s2_tokens_compact = [cell_id.to_token() for cell_id in covering.cell_ids()]
        s2_features = [] 
        for s2_token_compact in tqdm(s2_tokens_compact, desc="Compacting cells "):
            s2_id_compact = s2.CellId.from_token(s2_token_compact)
            s2_cell = s2.Cell(s2_id_compact)    
            # Get the vertices of the cell (4 vertices for a rectangular cell)
            vertices = [s2_cell.get_vertex(i) for i in range(4)]
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
            resolution = s2_id_compact.level()
            num_edges = 4
            s2_feature = geodesic_dggs_to_feature("s2",s2_token_compact,resolution,cell_polygon,num_edges)   
            s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
    
def s2compact_cli():
    """
    Command-line interface for s2compact.
    """
    parser = argparse.ArgumentParser(description="Compact S2 in a GeoJSON file containing a s2_token property named 's2'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input S2 in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = s2compact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "s2_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


def s2_expand(s2_ids, resolution):
    uncopmpacted_cells = []
    for s2_id in s2_ids:
        if s2_id.level() >= resolution:
            uncopmpacted_cells.append(s2_id)
        else:
            uncopmpacted_cells.extend(s2_id.children(resolution))  # Expand to the target level

    return uncopmpacted_cells

def s2expand(geojson_data,resolution):
    s2_tokens = [feature["properties"]["s2"] for feature in geojson_data.get("features", []) if "s2" in feature.get("properties", {})]
    s2_ids = [s2.CellId.from_token(token) for token in s2_tokens]
    if s2_ids:
        s2_ids_expand = s2_expand(s2_ids, resolution)
        s2_tokens_expand = [s2_id_expand.to_token() for s2_id_expand in s2_ids_expand]
        s2_features = [] 
        for s2_token_expand in tqdm(s2_tokens_expand, desc="Expanding cells "):
            s2_id_expand = s2.CellId.from_token(s2_token_expand)
            s2_cell_expand = s2.Cell(s2_id_expand)    
            # Get the vertices of the cell (4 vertices for a rectangular cell)
            vertices = [s2_cell_expand.get_vertex(i) for i in range(4)]
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
            resolution = s2_id_expand.level()
            num_edges = 4
            s2_feature = geodesic_dggs_to_feature("s2",s2_token_expand,resolution,cell_polygon,num_edges)   
            s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
    
def s2expand_cli():
    """
    Command-line interface for s2expand.
    """
    parser = argparse.ArgumentParser(description="expand S2 in a GeoJSON file containing a s2_token property named 's2'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of H3 to be expanded [0..30]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input S2 in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    resolution = args.resolution

    if resolution < 0 or resolution > 30:
        print(f"Please select a resolution in [0..30] range and try again ")
        return
    
    geojson_features = s2expand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"s2_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# rHEALPix
#################
def rhealpix_compact(rhealpix_dggs, rhealpix_ids):
    """
    Fully compacts rHEALPix cell IDs by replacing fully populated parent cells with their parent.
    Iterates until no further compaction is possible.
    """
    rhealpix_ids = set(rhealpix_ids)  # Remove duplicates
    
    # Main loop for compaction
    while True:
        grouped_rhealpix_ids = defaultdict(set)
        
        # Group cells by their parent
        for rhealpix_id in rhealpix_ids:
            if len(rhealpix_id) > 1:  # Ensure there's a valid parent
                parent = rhealpix_id[:-1]
                grouped_rhealpix_ids[parent].add(rhealpix_id)
        
        new_rhealpix_ids = set(rhealpix_ids)
        changed = False
        
        # Check if we can replace children with parent
        for parent, children in grouped_rhealpix_ids.items():
            parent_uids = (parent[0],) + tuple(map(int, parent[1:]))  # Assuming parent is a string like 'A0'
            parent_cell = rhealpix_dggs.cell(parent_uids)  # Retrieve the parent cell object
            
            # Generate the subcells for the parent at the next resolution
            subcells_at_next_res = set(str(subcell) for subcell in parent_cell.subcells())  # Collect subcells as strings
            
            # Check if the current children match the subcells at the next resolution
            if children == subcells_at_next_res:
                new_rhealpix_ids.difference_update(children)  # Remove children
                new_rhealpix_ids.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        rhealpix_ids = new_rhealpix_ids  # Continue compacting
    
    return sorted(rhealpix_ids)  # Sorted for consistency


def rhealpixcompact(rhealpix_dggs,geojson_data):
    rhealpix_ids = [feature["properties"]["rhealpix"] for feature in geojson_data.get("features", []) if "rhealpix" in feature.get("properties", {})]
    rhealpix_ids_compact = rhealpix_compact(rhealpix_dggs,rhealpix_ids)
    rhealpix_features = [] 
    for rhealpix_id_compact in tqdm(rhealpix_ids_compact, desc="Compacting cells "):  
        rhealpix_uids = (rhealpix_id_compact[0],) + tuple(map(int, rhealpix_id_compact[1:]))       
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        # if rhealpix_cell:
        resolution = rhealpix_cell.resolution        
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
        num_edges = 4
        if rhealpix_cell.ellipsoidal_shape() == 'dart':
            num_edges = 3
        rhealpix_feature = geodesic_dggs_to_feature("rhealpix",rhealpix_id_compact,resolution,cell_polygon,num_edges)   
        rhealpix_features.append(rhealpix_feature)

    return {
        "type": "FeatureCollection",
        "features": rhealpix_features
    }
           
def rhealpixcompact_cli():
    """
    Command-line interface for rhealpixcompact.
    """
    parser = argparse.ArgumentParser(description="Compact Rhealpix in a GeoJSON file containing a rhealpix ID property named 'rhealpix'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input Rhealpix in GeoJSON"
    )
    rhealpix_dggs = RHEALPixDGGS()

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = rhealpixcompact(rhealpix_dggs,geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "rhealpix_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
        

def rhealpix_expand(rhealpix_dggs, rhealpix_ids, resolution):
    expand_cells = []
    for rhealpix_id in rhealpix_ids: 
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))       
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        cell_resolution = rhealpix_cell.resolution  

        if cell_resolution >= resolution:
            expand_cells.append(rhealpix_cell)
        else:
            expand_cells.extend(rhealpix_cell.subcells(resolution))  # Expand to the target level
    return expand_cells


def rhealpixexpand(rhealpix_dggs,geojson_data,resolution):
    rhealpix_ids = [feature["properties"]["rhealpix"] for feature in geojson_data.get("features", []) if "rhealpix" in feature.get("properties", {})]
    rhealpix_cells_expand = rhealpix_expand(rhealpix_dggs,rhealpix_ids,resolution)
    rhealpix_features = [] 
    for rhealpix_cell_expand in rhealpix_cells_expand:    
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell_expand)
        rhealpix_id_expand = str(rhealpix_cell_expand)
        num_edges = 4
        if rhealpix_cell_expand.ellipsoidal_shape() == 'dart':
            num_edges = 3
        rhealpix_feature = geodesic_dggs_to_feature("rhealpix",rhealpix_id_expand,resolution,cell_polygon,num_edges)   
        rhealpix_features.append(rhealpix_feature)

    return {
        "type": "FeatureCollection",
        "features": rhealpix_features
    }
           
def rhealpixexpand_cli():
    """
    Command-line interface for rhealpixexpand.
    """
    parser = argparse.ArgumentParser(description="Uccompact Rhealpix in a GeoJSON file containing a rhealpix ID property named 'rhealpix'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input Rhealpix in GeoJSON"
    )
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of Rhealpix to be expanded [0..15]")
    
    rhealpix_dggs = RHEALPixDGGS()

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    
    if resolution < 0 or resolution > 15:
        print(f"Please select a resolution in [0..15] range and try again ")
        return


    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = rhealpixexpand(rhealpix_dggs,geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"rhealpix_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# ISEA4T
#################
def get_isea4t_cell_children(isea4t_dggs, isea4t_cell, resolution):
    if (platform.system() == 'Windows'): 
        """Recursively expands a DGGS cell until all children reach the desired resolution."""
        cell_id = isea4t_cell.get_cell_id()
        cell_resolution = len(cell_id) - 2  

        if cell_resolution >= resolution:
            return [isea4t_cell]  # Base case: return the cell if it meets/exceeds resolution

        expanded_cells = []
        children = isea4t_dggs.get_dggs_cell_children(isea4t_cell)

        for child in children:
            expanded_cells.extend(get_isea4t_cell_children(isea4t_dggs, child, resolution))

        return expanded_cells

def isea4t_compact(isea4t_dggs, isea4t_ids):
    if (platform.system() == 'Windows'):     
        isea4t_ids = set(isea4t_ids)  # Remove duplicates
        # Main loop for compaction
        while True:
            grouped_isea4t_ids = defaultdict(set)            
            # Group cells by their parent
            for isea4t_id in isea4t_ids:
                if len(isea4t_id) > 2:  # Ensure there's a valid parent
                    parent = isea4t_id[:-1]
                    grouped_isea4t_ids[parent].add(isea4t_id)
            
            new_isea4t_ids = set(isea4t_ids)
            changed = False
            
            # Check if we can replace children with parent
            for parent, children in grouped_isea4t_ids.items():
                parent_cell = DggsCell(parent)
                # Generate the subcells for the parent at the next resolution
                children_at_next_res = set(child.get_cell_id() for child in isea4t_dggs.get_dggs_cell_children(parent_cell))  # Collect subcells as strings
                
                # Check if the current children match the subcells at the next resolution
                if children == children_at_next_res:
                    new_isea4t_ids.difference_update(children)  # Remove children
                    new_isea4t_ids.add(parent)  # Add the parent
                    changed = True  # A change occurred
            
            if not changed:
                break  # Stop if no more compaction is possible
            isea4t_ids = new_isea4t_ids  # Continue compacting
        
        return sorted(isea4t_ids)  # Sorted for consistency

def isea4tcompact(isea4t_dggs,geojson_data):
    if (platform.system() == 'Windows'):
        isea4t_ids = [feature["properties"]["isea4t"] for feature in geojson_data.get("features", []) if "isea4t" in feature.get("properties", {})]
        isea4t_ids_compact = isea4t_compact(isea4t_dggs,isea4t_ids)
        isea4t_features = [] 
        for isea4t_id_compact in tqdm(isea4t_ids_compact, desc="Compacting cells "):    
            isea4t_cell_compact = DggsCell(isea4t_id_compact)
            cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(isea4t_cell_compact,ShapeStringFormat.WKT)
            cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
            if isea4t_id_compact.startswith('00') or isea4t_id_compact.startswith('09') or isea4t_id_compact.startswith('14')\
                or isea4t_id_compact.startswith('04') or isea4t_id_compact.startswith('19'):
                cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
            
            resolution = len(isea4t_id_compact) -2
            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
            num_edges = 3
            isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id_compact,resolution,cell_polygon,num_edges)   
            isea4t_features.append(isea4t_feature)

        return {
            "type": "FeatureCollection",
            "features": isea4t_features
        }

def isea4tcompact_cli():
    if (platform.system() == 'Windows'):  
        """
        Command-line interface for isea4tcompact.
        """
        isea4t_dggs = Eaggr(Model.ISEA4T)
        parser = argparse.ArgumentParser(description="Compact isea4t in a GeoJSON file containing a ISEA4T ID property named 'isea4t'")
        parser.add_argument(
            '-geojson', '--geojson', type=str, required=True, help="Input ISEA4T in GeoJSON"
        )

        args = parser.parse_args()
        geojson = args.geojson

        if not os.path.exists(geojson):
            print(f"Error: The file {geojson} does not exist.")
            return

        with open(geojson, 'r', encoding='utf-8') as f:
            geojson_data = json.load(f)
        
        geojson_features = isea4tcompact(isea4t_dggs,geojson_data)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = "isea4t_compacted.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")
        
        
def isea4t_expand(isea4t_dggs, isea4t_ids, resolution):
    """Expands a list of DGGS cells to the target resolution."""
    if platform.system() == 'Windows':       
        expand_cells = []
        for isea4t_id in isea4t_ids: 
            isea4t_cell = DggsCell(isea4t_id)
            expand_cells.extend(get_isea4t_cell_children(isea4t_dggs, isea4t_cell, resolution))
        return expand_cells

def isea4texpand(isea4t_dggs,geojson_data,resolution):
    if (platform.system() == 'Windows'):  
        isea4t_ids = [feature["properties"]["isea4t"] for feature in geojson_data.get("features", []) if "isea4t" in feature.get("properties", {})]
        isea4t_cells_expand = isea4t_expand(isea4t_dggs,isea4t_ids,resolution)
        isea4t_features = [] 
        for isea4t_cell_expand in tqdm(isea4t_cells_expand, desc="Expanding cells "):    
            cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(isea4t_cell_expand,ShapeStringFormat.WKT)
            cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
            isea4t_id_expand = isea4t_cell_expand.get_cell_id()
            if isea4t_id_expand.startswith('00') or isea4t_id_expand.startswith('09') or isea4t_id_expand.startswith('14')\
                or isea4t_id_expand.startswith('04') or isea4t_id_expand.startswith('19'):
                cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
                
            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
            num_edges = 3
            isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id_expand,resolution,cell_polygon,num_edges)   
            isea4t_features.append(isea4t_feature)

        return {
            "type": "FeatureCollection",
            "features": isea4t_features
        }
            
def isea4texpand_cli():
    if (platform.system() == 'Windows'):
        """
        Command-line interface for isea4texpand.
        """
        isea4t_dggs = Eaggr(Model.ISEA4T)
        parser = argparse.ArgumentParser(description="Uccompact OpenEaggr ISEA4T in a GeoJSON file containing an ISEA4T ID property named 'isea4t'")
        parser.add_argument(
            '-geojson', '--geojson', type=str, required=True, help="Input OpenEaggr ISEA4T in GeoJSON"
        )
        parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of OpenEaggr ISEA4T to be expanded [0..25]")


        args = parser.parse_args()
        geojson = args.geojson
        resolution = args.resolution
        
        # actual resolution range: [0..39]
        if resolution < 0 or resolution > 25:
            print(f"Please select a resolution in [0..25] range and try again ")
            return

        if not os.path.exists(geojson):
            print(f"Error: The file {geojson} does not exist.")
            return

        with open(geojson, 'r', encoding='utf-8') as f:
            geojson_data = json.load(f)
        
        geojson_features = isea4texpand(isea4t_dggs,geojson_data,resolution)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"isea4t_{resolution}_expanded.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")


#################
# ISEA3H
#################
def get_isea3h_cell_children(isea3h_dggs, isea3h_cell, resolution):
    if (platform.system() == 'Windows'): 
        """Recursively expands a DGGS cell until all children reach the desired resolution."""
        isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)      
        cell_accuracy = isea3h2point._accuracy            
        cell_resolution  = isea3h_accuracy_res_dict.get(cell_accuracy)
                    
        if cell_resolution >= resolution:
            return [isea3h_cell]  # Base case: return the cell if it meets/exceeds resolution

        expanded_cells = []
        children = isea3h_dggs.get_dggs_cell_children(isea3h_cell)

        for child in children:
            expanded_cells.extend(get_isea3h_cell_children(isea3h_dggs, child, resolution))

        return expanded_cells

def isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell):
    if (platform.system() == 'Windows'):
        cell_to_shape = isea3h_dggs.convert_dggs_cell_outline_to_shape_string(isea3h_cell,ShapeStringFormat.WKT)
        if cell_to_shape:
            coordinates_part = cell_to_shape.replace("POLYGON ((", "").replace("))", "")
            coordinates = []
            for coord_pair in coordinates_part.split(","):
                lon, lat = map(float, coord_pair.strip().split())
                coordinates.append([lon, lat])

            # Ensure the polygon is closed (first and last point must be the same)
            if coordinates[0] != coordinates[-1]:
                coordinates.append(coordinates[0])

        cell_polygon = Polygon(coordinates)
        fixed_polygon = fix_polygon(cell_polygon)    
        return fixed_polygon

def isea3h_compact(isea3h_dggs, isea3h_ids):
    isea3h_ids = set(isea3h_ids)  # Remove duplicates
    cell_cache = {cell_id: DggsCell(cell_id) for cell_id in isea3h_ids}

    while True:
        grouped_by_parent = defaultdict(set)

        # Group cells by *all* their parents
        for cell_id in isea3h_ids:
            cell = cell_cache[cell_id]
            try:
                parents = isea3h_dggs.get_dggs_cell_parents(cell)
            except Exception as e:
                print(f"Error getting parents for {cell_id}: {e}")
                continue

            for parent in parents:
                parent_id = parent.get_cell_id()
                grouped_by_parent[parent_id].add(cell_id)

        new_isea3h_ids = set(isea3h_ids)
        changed = False

        for parent_id, children_ids in grouped_by_parent.items():
            parent_cell = DggsCell(parent_id)
            try:
                expected_children = set(
                    child.get_cell_id()
                    for child in isea3h_dggs.get_dggs_cell_children(parent_cell)
                )
            except Exception as e:
                print(f"Error getting children for parent {parent_id}: {e}")
                continue

            # Check for full match: only then compact
            if children_ids == expected_children:
                new_isea3h_ids.difference_update(children_ids)
                new_isea3h_ids.add(parent_id)
                cell_cache[parent_id] = parent_cell
                changed = True
            ########## 
            else:
                # Keep original children if they don't fully match expected subcells
                new_isea3h_ids.update(children_ids)


        if not changed:
            break  # Fully compacted

        isea3h_ids = new_isea3h_ids

    return sorted(isea3h_ids)

def isea3hcompact(isea3h_dggs,geojson_data):
    if (platform.system() == 'Windows'):
        isea3h_ids = [feature["properties"]["isea3h"] for feature in geojson_data.get("features", []) if "isea3h" in feature.get("properties", {})]
        isea3h_ids_compact = isea3h_compact(isea3h_dggs,isea3h_ids)
        isea3h_features = [] 
        for isea3h_id_compact in tqdm(isea3h_ids_compact, desc="Compacting cells "):    
            isea3h_cell = DggsCell(isea3h_id_compact)
            
            cell_polygon = isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell)
            isea3h_id = isea3h_cell.get_cell_id()
            cell_centroid = cell_polygon.centroid
            center_lat =  round(cell_centroid.y, 7)
            center_lon = round(cell_centroid.x, 7)
            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            
            isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)      
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
            isea3h_features.append(isea3h_feature)
    
    return {
        "type": "FeatureCollection",
        "features": isea3h_features,
    }

def isea3hcompact_cli():
    if (platform.system() == 'Windows'):  
        isea3h_dggs = Eaggr(Model.ISEA3H)
        """
        Command-line interface for isea3hcompact.
        """
        parser = argparse.ArgumentParser(description="Compact ISEA3H in a GeoJSON file containing a ISEA3H ID property named 'isea3h'")
        parser.add_argument(
            '-geojson', '--geojson', type=str, required=True, help="Input ISEA3H in GeoJSON"
        )

        args = parser.parse_args()
        geojson = args.geojson

        if not os.path.exists(geojson):
            print(f"Error: The file {geojson} does not exist.")
            return

        with open(geojson, 'r', encoding='utf-8') as f:
            geojson_data = json.load(f)
        
        geojson_features = isea3hcompact(isea3h_dggs,geojson_data)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = "isea3h_compacted.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")
        
        
def isea3h_expand(isea3h_dggs, isea3h_ids, resolution):
    """Expands a list of DGGS cells to the target resolution."""
    if platform.system() == 'Windows':       
        expand_cells = []
        for isea3h_id in isea3h_ids: 
            isea3h_cell = DggsCell(isea3h_id)
            expand_cells.extend(get_isea3h_cell_children(isea3h_dggs, isea3h_cell, resolution))
        return expand_cells

def isea3hexpand(isea3h_dggs,geojson_data,resolution):
    if (platform.system() == 'Windows'):  
        isea3h_ids = [feature["properties"]["isea3h"] for feature in geojson_data.get("features", []) if "isea3h" in feature.get("properties", {})]
        isea3h_cells_expand = isea3h_expand(isea3h_dggs,isea3h_ids,resolution)
        isea3h_features = [] 
        for isea3h_cell_expand in tqdm(isea3h_cells_expand, desc="Expanding cells "):    
            cell_polygon = isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell_expand)
           
            isea3h_id = isea3h_cell_expand.get_cell_id()
            cell_centroid = cell_polygon.centroid
            center_lat =  round(cell_centroid.y, 7)
            center_lon = round(cell_centroid.x, 7)
            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            
            isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell_expand)      
            cell_accuracy = isea3h2point._accuracy
                
            avg_edge_len = cell_perimeter / 6
            cell_resolution  = isea3h_accuracy_res_dict.get(cell_accuracy)
            
            if (cell_resolution == 0): # icosahedron faces at resolution = 0
                avg_edge_len = cell_perimeter / 3
            
            if cell_accuracy == 0.0:
                if round(avg_edge_len,2) == 0.06:
                    cell_resolution = 33
                elif round(avg_edge_len,2) == 0.03:
                    cell_resolution = 34
                elif round(avg_edge_len,2) == 0.02:
                    cell_resolution = 35
                elif round(avg_edge_len,2) == 0.01:
                    cell_resolution = 36
                
                elif round(avg_edge_len,3) == 0.007:
                    cell_resolution = 37
                elif round(avg_edge_len,3) == 0.004:
                    cell_resolution = 38
                elif round(avg_edge_len,3) == 0.002:
                    cell_resolution = 39
                elif round(avg_edge_len,3) <= 0.001:
                    cell_resolution = 40
        
            isea3h_feature = {
                "type": "Feature",
                "geometry": mapping(cell_polygon),
                "properties": {
                        "isea3h": isea3h_id,
                        "resolution": cell_resolution,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "avg_edge_len": round(avg_edge_len,3),
                        "cell_area": cell_area
                        }
            }
            isea3h_features.append(isea3h_feature)
    
    return {
        "type": "FeatureCollection",
        "features": isea3h_features,
    }

            
def isea3hexpand_cli():
    if (platform.system() == 'Windows'):
        """
        Command-line interface for isea3hexpand.
        """
        isea3h_dggs = Eaggr(Model.ISEA3H)
        parser = argparse.ArgumentParser(description="Uccompact OpenEaggr ISEA3H in a GeoJSON file containing an ISEA3H ID property named 'isea3h'")
        parser.add_argument(
            '-geojson', '--geojson', type=str, required=True, help="Input OpenEaggr ISEA3H in GeoJSON"
        )
        parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of OpenEaggr ISEA3H to be expanded [0..32]")


        args = parser.parse_args()
        geojson = args.geojson
        resolution = args.resolution
        
        # actual resolution range: [0..40]
        if resolution < 0 or resolution > 32:
            print(f"Please select a resolution in [0..32] range and try again ")
            return

        if not os.path.exists(geojson):
            print(f"Error: The file {geojson} does not exist.")
            return

        with open(geojson, 'r', encoding='utf-8') as f:
            geojson_data = json.load(f)
        
        geojson_features = isea3hexpand(isea3h_dggs,geojson_data,resolution)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"isea3h_{resolution}_expanded.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")



#################
# EASE
#################
def ease_compact (ease_ids):   
    ease_ids = set(ease_ids)  # Remove duplicates

    while True:
        grouped_ease_ids = defaultdict(set)

        # Group cells by their parent
        for ease_id in ease_ids:
            match = re.match(r"L(\d+)\.(.+)", ease_id)  # Extract resolution level & ID
            if not match:
                continue  # Skip invalid IDs
            
            resolution = int(match.group(1))
            base_id = match.group(2)

            if resolution == 0:
                continue  # L0 has no parent

            # Determine the parent by removing the last section
            parent = f"L{resolution-1}." + ".".join(base_id.split(".")[:-1])
            # print (f"parent: {parent}")
            grouped_ease_ids[parent].add(ease_id)

        new_ease_ids = set(ease_ids)
        changed = False

        # Check if we can replace children with their parent
        for parent, children in grouped_ease_ids.items():
            # print (f"children: {children}")
            match = re.match(r"L(\d+)\..+", parent)
            if not match:
                continue  # Skip invalid parents

            resolution = int(match.group(1))
            children_at_next_res = set(_parent_to_children(parent, resolution+1))  # Ensure correct format
            # If all expected children are present, replace them with the parent
            if children == children_at_next_res:  
                new_ease_ids.difference_update(children)
                new_ease_ids.add(parent)
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        ease_ids = new_ease_ids  # Continue compacting
    
    return sorted(ease_ids)  # Sorted for consistency

def easecompact(geojson_data):
    ease_ids = [feature["properties"]["ease"] for feature in geojson_data.get("features", []) if "ease" in feature.get("properties", {})]
    ease_cells_compact = ease_compact(ease_ids)
    ease_features = [] 
    for ease_cell_compact in tqdm(ease_cells_compact, desc="Compacting cells "):    
        level = int(ease_cell_compact[1])  # Get the level (e.g., 'L0' -> 0)
        # Get level specs
        level_spec = levels_specs[level]
        n_row = level_spec["n_row"]
        n_col = level_spec["n_col"]
            
        geo = grid_ids_to_geos([ease_cell_compact])
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
            ease_feature = geodesic_dggs_to_feature("ease",ease_cell_compact,resolution,cell_polygon,num_edges)   
            ease_features.append(ease_feature)

    return {
        "type": "FeatureCollection",
        "features": ease_features
    }


def easecompact_cli():
    """
    Command-line interface for easecompact.
    """
    parser = argparse.ArgumentParser(description="Compact EASE in a GeoJSON file containing a EASE ID property named 'ease'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input EASE in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = easecompact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "ease_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}") 

def ease_expand(ease_ids, resolution):
    uncopmpacted_cells = []
    for ease_id in ease_ids:
        ease_resolution = int(ease_id[1])
        if ease_resolution >= resolution:
            uncopmpacted_cells.append(ease_id)
        else:
            uncopmpacted_cells.extend(_parent_to_children(ease_id, ease_resolution+1))  # Expand to the target level

    return uncopmpacted_cells


def easeexpand(geojson_data,resolution):
    ease_ids = [feature["properties"]["ease"] for feature in geojson_data.get("features", []) if "ease" in feature.get("properties", {})]
    ease_cells_expand = ease_expand(ease_ids,resolution)
    ease_features = [] 
    for ease_cell_expand in tqdm(ease_cells_expand, desc="Expanding cells "):    
        level = int(ease_cell_expand[1])  # Get the level (e.g., 'L0' -> 0)
        # Get level specs
        level_spec = levels_specs[level]
        n_row = level_spec["n_row"]
        n_col = level_spec["n_col"]
            
        geo = grid_ids_to_geos([ease_cell_expand])
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
            ease_feature = geodesic_dggs_to_feature("ease",ease_cell_expand,resolution,cell_polygon,num_edges)   
            ease_features.append(ease_feature)

    return {
        "type": "FeatureCollection",
        "features": ease_features
    }

def easeexpand_cli():
    """
    Command-line interface for easeexpand.
    """
    parser = argparse.ArgumentParser(description="expand EASE in a GeoJSON file containing a EASE ID property named 'ease'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input EASE in GeoJSON"
    )

    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of EASE to be expanded [0..6]")


    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution < 0 or resolution > 6:
        print(f"Please select a resolution in [0..6] range and try again ")
        return


    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = easeexpand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "ease_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}") 


#################
# QTM
#################
def qtm_compact(qtm_ids):
    qtm_ids = set(qtm_ids)  # Remove duplicates
    # Main loop for compaction
    while True:
        grouped_qtm_ids = defaultdict(set)        
        # Group cells by their parent
        for qtm_id in qtm_ids:
            parent = qtm.qtm_parent(qtm_id)
            grouped_qtm_ids[parent].add(qtm_id)

        new_qtm_ids = set(qtm_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_qtm_ids.items():
            next_resolution = len(parent)+1          
            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(childcell for childcell in qtm.qtm_children(parent,next_resolution)) 
            
            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_qtm_ids.difference_update(children)  # Remove children
                new_qtm_ids.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        qtm_ids = new_qtm_ids  # Continue compacting
    
    return sorted(qtm_ids)  # Sorted for consistency

def qtmcompact(geojson_data):
    qtm_ids = [feature["properties"]["qtm"] for feature in geojson_data.get("features", []) if "qtm" in feature.get("properties", {})]
    qtm_ids_compact = qtm_compact(qtm_ids)
    qtm_features = [] 
    for qtm_id_compact in tqdm(qtm_ids_compact, desc="Compacting cells "):  
        facet = qtm.qtm_id_to_facet(qtm_id_compact)
        cell_polygon = qtm.constructGeometry(facet)    
        cell_resolution = len(qtm_id_compact)
        num_edges = 3

        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id_compact,cell_resolution,cell_polygon,num_edges)  
        qtm_features.append(qtm_feature)

    return {
        "type": "FeatureCollection",
        "features": qtm_features
    }
    
def qtmcompact_cli():
    """
    Command-line interface for qtmcompact.
    """
    parser = argparse.ArgumentParser(description="Compact QTM in GeoJSON format containing a property named 'qtm'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input QTM in GeoJSON format"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = qtmcompact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "qtm_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


def qtm_expand(qtm_ids, resolution):
    expand_cells = []
    for qtm_id in qtm_ids: 
        cell_resolution = len(qtm_id)
        if cell_resolution >= resolution:
            expand_cells.append(qtm_id)
        else:
            expand_cells.extend(qtm.qtm_children(qtm_id,resolution))  # Expand to the target level
    return expand_cells
      
def qtmexpand(geojson_data,resolution):
    qtm_ids = [feature["properties"]["qtm"] for feature in geojson_data.get("features", []) if "qtm" in feature.get("properties", {})]
    qtm_ids_expand = qtm_expand(qtm_ids, resolution)
    qtm_features = [] 
    for qtm_id_expand in tqdm(qtm_ids_expand, desc="Expanding cells "):
        facet = qtm.qtm_id_to_facet(qtm_id_expand)
        cell_polygon = qtm.constructGeometry(facet)    
        cell_resolution = len(qtm_id_expand)
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id_expand,cell_resolution,cell_polygon,num_edges)  
        qtm_features.append(qtm_feature)

    return {
        "type": "FeatureCollection",
        "features": qtm_features
    }
    
    
def qtmexpand_cli():
    """
    Command-line interface for qtmexpand.
    """
    parser = argparse.ArgumentParser(description="expand QTM in GeoJSON format containing a property named 'qtm'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of qtm to be expanded [1..24]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input qtm in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution < 1 or resolution > 24:
        print(f"Please select a resolution in [1..24] and try again ")
        return

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = qtmexpand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"qtm_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# OLC
#################
def olc_compact(olc_ids):
    olc_ids = set(olc_ids)  # Remove duplicates
    # Main loop for compaction
    while True:
        grouped_olc_ids = defaultdict(set)        
        # Group cells by their parent
        for olc_id in olc_ids:
            parent = olc.olc_parent(olc_id)
            grouped_olc_ids[parent].add(olc_id)

        new_olc_ids = set(olc_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_olc_ids.items():
            coord = olc.decode(parent) 
            coord_len = coord.codeLength 
            if coord_len <= 10:
                next_resolution = coord_len + 2
            else:
                next_resolution = coord_len + 1
          
            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(childcell for childcell in olc.olc_children(parent,next_resolution)) 
            
            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_olc_ids.difference_update(children)  # Remove children
                new_olc_ids.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        olc_ids = new_olc_ids  # Continue compacting
    
    return sorted(olc_ids)  # Sorted for consistency

def olccompact(geojson_data):
    olc_ids = [feature["properties"]["olc"] for feature in geojson_data.get("features", []) if "olc" in feature.get("properties", {})]
    olc_ids_compact = olc_compact(olc_ids)
    olc_features = [] 
    for olc_id_compact in tqdm(olc_ids_compact, desc="Compacting cells "):  
        coord = olc.decode(olc_id_compact)    
        if coord:
            # Create the bounding box coordinates for the polygon
            min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
            max_lat, max_lon = coord.latitudeHi, coord.longitudeHi        
            cell_resolution = coord.codeLength 

            # Define the polygon based on the bounding box
            cell_polygon = Polygon([
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat]   # Closing the polygon (same as the first point)
            ])
            olc_feature = graticule_dggs_to_feature("olc",olc_id_compact,cell_resolution,cell_polygon)   
            olc_features.append(olc_feature)

    return {
        "type": "FeatureCollection",
        "features": olc_features
    }
    
def olccompact_cli():
    """
    Command-line interface for olccompact.
    """
    parser = argparse.ArgumentParser(description="Compact olc in a GeoJSON file containing a property named 'olc'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input olc in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = olccompact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "olc_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


def olc_expand(olc_ids, resolution):
    expand_cells = []
    for olc_id in olc_ids: 
        coord = olc.decode(olc_id)  
        cell_resolution = coord.codeLength
        if cell_resolution >= resolution:
            expand_cells.append(olc_id)
        else:
            expand_cells.extend(olc.olc_children(olc_id,resolution))  # Expand to the target level
    return expand_cells

def olcexpand(geojson_data,resolution):
    olc_ids = [feature["properties"]["olc"] for feature in geojson_data.get("features", []) if "olc" in feature.get("properties", {})]
    olc_ids_expand = olc_expand(olc_ids, resolution)
    olc_features = [] 
    for olc_id_expand in tqdm(olc_ids_expand, desc="Expanding cells "):
        coord = olc.decode(olc_id_expand)    
        if coord:
            # Create the bounding box coordinates for the polygon
            min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
            max_lat, max_lon = coord.latitudeHi, coord.longitudeHi        

            # Define the polygon based on the bounding box
            cell_polygon = Polygon([
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat]   # Closing the polygon (same as the first point)
            ])
            olc_feature = graticule_dggs_to_feature("olc",olc_id_expand,resolution,cell_polygon)   
            olc_features.append(olc_feature)

    return {
        "type": "FeatureCollection",
        "features": olc_features
    }
    
    
def olcexpand_cli():
    """
    Command-line interface for olcexpand.
    """
    parser = argparse.ArgumentParser(description="expand olc in a GeoJSON file containing a property named 'olc'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of olc to be expanded [2, 4, 6, 8, 10..15]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input olc in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution not in [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]:
        print(f"Please select a resolution in [2, 4, 6, 8, 10..15] and try again ")
        return

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = olcexpand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"olc_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# Geohash
#################
def geohash_compact(geohash_ids):
    geohash_ids = set(geohash_ids)  # Remove duplicates
    # Main loop for compaction
    while True:
        grouped_geohash_ids = defaultdict(set)        
        # Group cells by their parent
        for geohash_id in geohash_ids:
            parent = geohash.geohash_parent(geohash_id)
            grouped_geohash_ids[parent].add(geohash_id)

        new_geohash_ids = set(geohash_ids)
        changed = False


        # Check if we can replace children with parent
        for parent, children in grouped_geohash_ids.items():
            parent_resolution =  len(parent)
            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(childcell for childcell in geohash.geohash_children(parent,parent_resolution+1)) 
            
            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_geohash_ids.difference_update(children)  # Remove children
                new_geohash_ids.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        geohash_ids = new_geohash_ids  # Continue compacting
    
    return sorted(geohash_ids)  # Sorted for consistency

def geohashcompact(geojson_data):
    geohash_ids = [feature["properties"]["geohash"] for feature in geojson_data.get("features", []) if "geohash" in feature.get("properties", {})]
    geohash_ids_compact = geohash_compact(geohash_ids)
    geohash_features = [] 
    for geohash_id_compact in tqdm(geohash_ids_compact, desc="Compacting cells "):  
        cell_polygon = geohash_to_polygon(geohash_id_compact)
        cell_resolution =  len(geohash_id_compact)
        geohash_feature = graticule_dggs_to_feature("geohash",geohash_id_compact,cell_resolution,cell_polygon)   
        geohash_features.append(geohash_feature)

    return {
        "type": "FeatureCollection",
        "features": geohash_features
    }
    
def geohashcompact_cli():
    """
    Command-line interface for geohashcompact.
    """
    parser = argparse.ArgumentParser(description="Compact Geohash in a GeoJSON file containing a property named 'geohash'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input Geohash in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = geohashcompact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "geohash_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


def geohash_expand(geohash_ids, resolution):
    expand_cells = []
    for geohash_id in geohash_ids: 
        cell_resolution = len(geohash_id)
        if cell_resolution >= resolution:
            expand_cells.append(geohash_id)
        else:
            expand_cells.extend(geohash.geohash_children(geohash_id,resolution))  # Expand to the target level
    return expand_cells

def geohashexpand(geojson_data,resolution):
    geohash_ids = [feature["properties"]["geohash"] for feature in geojson_data.get("features", []) if "geohash" in feature.get("properties", {})]
    geohash_ids_expand = geohash_expand(geohash_ids, resolution)
    geohash_features = [] 
    for geohash_id_expand in tqdm(geohash_ids_expand, desc="Expanding cells "):
        cell_polygon = geohash_to_polygon(geohash_id_expand)
        cell_resolution =  len(geohash_id_expand)
        geohash_feature = graticule_dggs_to_feature("geohash",geohash_id_expand,cell_resolution,cell_polygon)   
        geohash_features.append(geohash_feature)

    return {
        "type": "FeatureCollection",
        "features": geohash_features
    }
    
    
def geohashexpand_cli():
    """
    Command-line interface for geohashexpand.
    """
    parser = argparse.ArgumentParser(description="expand Geohash in a GeoJSON file containing a property named 'geohash'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of geohash to be expanded [1..10]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input geohash in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution < 1 or resolution > 10:
        print(f"Please select a resolution in [0..29] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = geohashexpand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"geohash_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# Tilecode
#################
def tilecode_compact(tilecode_ids):
    tilecode_ids = set(tilecode_ids)  # Remove duplicates
    
    # Main loop for compaction
    while True:
        grouped_tilecode_ids = defaultdict(set)
        
        # Group cells by their parent
        for tilecode_id in tilecode_ids:
            match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)    
            if match:  # Ensure there's a valid parent
                parent = tilecode.tilecode_parent(tilecode_id)
                grouped_tilecode_ids[parent].add(tilecode_id)

        new_tilecode_ids = set(tilecode_ids)
        changed = False
        
        # Check if we can replace children with parent
        for parent, children in grouped_tilecode_ids.items():
            # Generate the subcells for the parent at the next resolution
            match = re.match(r'z(\d+)x(\d+)y(\d+)', parent)    
            parent_resolution = int(match.group(1))

            childcells_at_next_res = set(childcell for childcell in tilecode.tilecode_children(parent,parent_resolution+1))  # Collect subcells as strings
            
            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_tilecode_ids.difference_update(children)  # Remove children
                new_tilecode_ids.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        tilecode_ids = new_tilecode_ids  # Continue compacting
    
    return sorted(tilecode_ids)  # Sorted for consistency

def tilecodecompact(geojson_data):
    tilecode_ids = [feature["properties"]["tilecode"] for feature in geojson_data.get("features", []) if "tilecode" in feature.get("properties", {})]
    tilecode_ids_compact = tilecode_compact(tilecode_ids)
    tilecode_features = [] 
    for tilecode_id_compact in tqdm(tilecode_ids_compact, desc="Compacting cells "):  
        match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id_compact)
        if not match:
            raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

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
            
            resolution = z
            tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id_compact,resolution,cell_polygon)   
            tilecode_features.append(tilecode_feature)

    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }
    
def tilecodecompact_cli():
    """
    Command-line interface for tilecodecompact.
    """
    parser = argparse.ArgumentParser(description="Compact Tilecode in a GeoJSON file containing a property named 'tilecode'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input Tilecode in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = tilecodecompact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "tilecode_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


def tilecode_expand(tilecode_ids, resolution):
    expand_cells = []
    for tilecode_id in tilecode_ids: 
        match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)
        if not match:
            raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")
        cell_resolution= int(match.group(1))
               
        if cell_resolution >= resolution:
            expand_cells.append(tilecode_id)
        else:
            expand_cells.extend(tilecode.tilecode_children(tilecode_id,resolution))  # Expand to the target level
    return expand_cells

def tilecodeexpand(geojson_data,resolution):
    tilecode_ids = [feature["properties"]["tilecode"] for feature in geojson_data.get("features", []) if "tilecode" in feature.get("properties", {})]
    tilecode_ids_expand = tilecode_expand(tilecode_ids, resolution)
    tilecode_features = [] 
    for tilecode_id_expand in tqdm(tilecode_ids_expand, desc="Expanding cells "):
        match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id_expand)
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
            
            resolution = z
            tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id_expand,resolution,cell_polygon)   
            tilecode_features.append(tilecode_feature)

    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }
    
    
def tilecodeexpand_cli():
    """
    Command-line interface for tilecodeexpand.
    """
    parser = argparse.ArgumentParser(description="expand Tilecode in a GeoJSON file containing a property named 'tilecode'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of Tilecode to be expanded [0..29]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input Tilecode in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution < 0 or resolution > 29:
        print(f"Please select a resolution in [0..29] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = tilecodeexpand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"tilecode_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# Quadkey
#################
def quadkey_compact(quadkey_ids):
    quadkey_ids = set(quadkey_ids)  # Remove duplicates
    
    # Main loop for compaction
    while True:
        grouped_quadkey_ids = defaultdict(set)
        
        # Group cells by their parent
        for quadkey_id in quadkey_ids:
            parent = tilecode.quadkey_parent(quadkey_id)
            grouped_quadkey_ids[parent].add(quadkey_id)

        new_quadkey_ids = set(quadkey_ids)
        changed = False


        # Check if we can replace children with parent
        for parent, children in grouped_quadkey_ids.items():
            parent_resolution = mercantile.quadkey_to_tile(parent).z

            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(childcell for childcell in tilecode.quadkey_children(parent,parent_resolution+1)) 
            
            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_quadkey_ids.difference_update(children)  # Remove children
                new_quadkey_ids.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        quadkey_ids = new_quadkey_ids  # Continue compacting
    
    return sorted(quadkey_ids)  # Sorted for consistency

def quadkeycompact(geojson_data):
    quadkey_ids = [feature["properties"]["quadkey"] for feature in geojson_data.get("features", []) if "quadkey" in feature.get("properties", {})]
    quadkey_ids_compact = quadkey_compact(quadkey_ids)
    quadkey_features = [] 
    for quadkey_id_compact in tqdm(quadkey_ids_compact, desc="Compacting cells "):  
        quadkey_id_compact_tile = mercantile.quadkey_to_tile(quadkey_id_compact)
        # Convert matched groups to integers
        z = quadkey_id_compact_tile.z
        x = quadkey_id_compact_tile.x
        y = quadkey_id_compact_tile.y

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
            
            cell_resolution = z
            quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id_compact,cell_resolution,cell_polygon)   
            quadkey_features.append(quadkey_feature)

    return {
        "type": "FeatureCollection",
        "features": quadkey_features
    }
    
def quadkeycompact_cli():
    """
    Command-line interface for quadkeycompact.
    """
    parser = argparse.ArgumentParser(description="Compact quadkey in a GeoJSON file containing a property named 'quadkey'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input quadkey in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = quadkeycompact(geojson_data)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = "quadkey_compacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


def quadkey_expand(quadkey_ids, resolution):
    expand_cells = []
    for quadkey_id in quadkey_ids: 
        quadkey_tile = mercantile.quadkey_to_tile(quadkey_id)
        cell_resolution = quadkey_tile.z
        if cell_resolution >= resolution:
            expand_cells.append(quadkey_id)
        else:
            expand_cells.extend(tilecode.quadkey_children(quadkey_id,resolution))  # Expand to the target level
    return expand_cells

def quadkeyexpand(geojson_data,resolution):
    quadkey_ids = [feature["properties"]["quadkey"] for feature in geojson_data.get("features", []) if "quadkey" in feature.get("properties", {})]
    quadkey_ids_expand = quadkey_expand(quadkey_ids, resolution)
    quadkey_features = [] 
    for quadkey_id_expand in tqdm(quadkey_ids_expand, desc="Expanding cells "):
        quadkey_id_expand_tile = mercantile.quadkey_to_tile(quadkey_id_expand)
        z = quadkey_id_expand_tile.z
        x = quadkey_id_expand_tile.x
        y = quadkey_id_expand_tile.y
        
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
            quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id_expand,resolution,cell_polygon)   
            quadkey_features.append(quadkey_feature)

    return {
        "type": "FeatureCollection",
        "features": quadkey_features
    }
    
    
def quadkeyexpand_cli():
    """
    Command-line interface for quadkeyexpand.
    """
    parser = argparse.ArgumentParser(description="expand quadkey in a GeoJSON file containing a property named 'quadkey'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of quadkey to be expanded [0..29]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input quadkey in GeoJSON"
    )

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution

    if resolution < 0 or resolution > 29:
        print(f"Please select a resolution in [0..29] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = quadkeyexpand(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"quadkey_{resolution}_expanded.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")

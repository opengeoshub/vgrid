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
from vgrid.conversion.dggs2geojson import h32geojson
from pyproj import Geod
geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID
from collections import defaultdict

from tqdm import tqdm
rhealpix_dggs = RHEALPixDGGS()
from vgrid.conversion.dggs2geojson import rhealpix_cell_to_polygon
#################
# H3
#################
def h3compact(geojson_data):
    h3_ids = list(set(feature["properties"]["h3"] for feature in geojson_data.get("features", []) if "h3" in feature.get("properties", {})))
    h3_ids_compact = h3.compact_cells(h3_ids)
    h3_features = [] 
    for h3_id_compact in tqdm(h3_ids_compact, desc="Processing cells "):  
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

def h3uncompact(geojson_data,resolution):
    h3_ids = [feature["properties"]["h3"] for feature in geojson_data.get("features", []) if "h3" in feature.get("properties", {})]
    h3_ids_uncompact = h3.uncompact_cells(h3_ids, resolution)
    h3_features = [] 
    for h3_id_uncompact in tqdm(h3_ids_uncompact, desc="Processing cells "):
        cell_boundary = h3.cell_to_boundary(h3_id_uncompact)   
        if cell_boundary:
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            num_edges = 6
            if (h3.is_pentagon(h3_id_uncompact)):
                num_edges = 5
            h3_feature = geodesic_dggs_to_feature("h3",h3_id_uncompact,resolution,cell_polygon,num_edges)   
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features
    }
    
def h3uncompact_cli():
    """
    Command-line interface for h3uncompact.
    """
    parser = argparse.ArgumentParser(description="Uncompact H3 in a GeoJSON file containing a property named 'h3'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of H3 to be uncompacted [0..15]")
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
    
    geojson_features = h3uncompact(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"h3_{resolution}_uncompacted.geojson"
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
        for s2_token_compact in tqdm(s2_tokens_compact, desc="Processing cells "):
            s2_id_compact = s2.CellId.from_token(s2_token_compact)
            s2_compact = s2.Cell(s2_id_compact)    
            # Get the vertices of the cell (4 vertices for a rectangular cell)
            vertices = [s2_compact.get_vertex(i) for i in range(4)]
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


def s2_uncompact(s2_ids, resolution):
    uncopmpacted_cells = []
    for s2_id in s2_ids:
        if s2_id.level() >= resolution:
            uncopmpacted_cells.append(s2_id)
        else:
            uncopmpacted_cells.extend(s2_id.children(resolution))  # Expand to the target level

    return uncopmpacted_cells

def s2uncompact(geojson_data,resolution):
    s2_tokens = [feature["properties"]["s2"] for feature in geojson_data.get("features", []) if "s2" in feature.get("properties", {})]
    s2_ids = [s2.CellId.from_token(token) for token in s2_tokens]
    if s2_ids:
        s2_ids_uncompact = s2_uncompact(s2_ids, resolution)
        s2_tokens_uncompact = [s2_id_uncompact.to_token() for s2_id_uncompact in s2_ids_uncompact]
        s2_features = [] 
        for s2_token_uncompact in tqdm(s2_tokens_uncompact, desc="Processing cells "):
            s2_id_uncompact = s2.CellId.from_token(s2_token_uncompact)
            s2_cell_uncompact = s2.Cell(s2_id_uncompact)    
            # Get the vertices of the cell (4 vertices for a rectangular cell)
            vertices = [s2_cell_uncompact.get_vertex(i) for i in range(4)]
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
            resolution = s2_id_uncompact.level()
            num_edges = 4
            s2_feature = geodesic_dggs_to_feature("s2",s2_token_uncompact,resolution,cell_polygon,num_edges)   
            s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
    
def s2uncompact_cli():
    """
    Command-line interface for s2uncompact.
    """
    parser = argparse.ArgumentParser(description="Uncompact S2 in a GeoJSON file containing a s2_token property named 's2'")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of H3 to be uncompacted [0..30]")
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
    
    geojson_features = s2uncompact(geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"s2_{resolution}_uncompacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


#################
# Rhealpix
#################
def rhealpix_compact(rhealpix_dggs, cells: list[str]) -> list[str]:
    """
    Fully compacts RHEALPix cell IDs by replacing fully populated parent cells with their parent.
    Iterates until no further compaction is possible.
    """
    cells = set(cells)  # Remove duplicates
    
    # Main loop for compaction
    while True:
        grouped_cells = defaultdict(set)
        
        # Group cells by their parent
        for cell in cells:
            if len(cell) > 1:  # Ensure there's a valid parent
                parent = cell[:-1]
                grouped_cells[parent].add(cell)
        
        new_cells = set(cells)
        changed = False
        
        # Check if we can replace children with parent
        for parent, children in grouped_cells.items():
            parent_uids = (parent[0],) + tuple(map(int, parent[1:]))  # Assuming parent is a string like 'A0'
            parent_cell = rhealpix_dggs.cell(parent_uids)  # Retrieve the parent cell object
            
            # Generate the subcells for the parent at the next resolution
            subcells_at_next_res = set(str(subcell) for subcell in parent_cell.subcells())  # Collect subcells as strings
            
            # Check if the current children match the subcells at the next resolution
            if children == subcells_at_next_res:
                new_cells.difference_update(children)  # Remove children
                new_cells.add(parent)  # Add the parent
                changed = True  # A change occurred
        
        if not changed:
            break  # Stop if no more compaction is possible
        cells = new_cells  # Continue compacting
    
    return sorted(cells)  # Sorted for consistency


def rhealpixcompact(rhealpix_dggs,geojson_data):
    rhealpix_ids = [feature["properties"]["rhealpix"] for feature in geojson_data.get("features", []) if "rhealpix" in feature.get("properties", {})]
    rhealpix_ids_compact = rhealpix_compact(rhealpix_dggs,rhealpix_ids)
    rhealpix_features = [] 
    for rhealpix_id_compact in tqdm(rhealpix_ids_compact, desc="Processing cells "):  
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
        

def rhealpix_uncompact(rhealpix_dggs, rhealpix_ids, resolution):
    uncompact_cells = []
    for rhealpix_id in rhealpix_ids:
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))       
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        cell_resolution = rhealpix_cell.resolution  

        if cell_resolution >= resolution:
            uncompact_cells.append(rhealpix_cell)
        else:
            uncompact_cells.extend(rhealpix_cell.subcells(resolution))  # Expand to the target level
    return uncompact_cells


def rhealpixuncompact(rhealpix_dggs,geojson_data,resolution):
    rhealpix_ids = [feature["properties"]["rhealpix"] for feature in geojson_data.get("features", []) if "rhealpix" in feature.get("properties", {})]
    rhealpix_cells_uncompact = rhealpix_uncompact(rhealpix_dggs,rhealpix_ids,resolution)
    rhealpix_features = [] 
    for rhealpix_cell_uncompact in tqdm(rhealpix_cells_uncompact, desc="Processing cells "):    
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell_uncompact)
        rhealpix_id_uncompact = str(rhealpix_cell_uncompact)
        num_edges = 4
        if rhealpix_cell_uncompact.ellipsoidal_shape() == 'dart':
            num_edges = 3
        rhealpix_feature = geodesic_dggs_to_feature("rhealpix",rhealpix_id_uncompact,resolution,cell_polygon,num_edges)   
        rhealpix_features.append(rhealpix_feature)

    return {
        "type": "FeatureCollection",
        "features": rhealpix_features
    }
           
def rhealpixuncompact_cli():
    """
    Command-line interface for rhealpixcompact.
    """
    parser = argparse.ArgumentParser(description="Uccompact Rhealpix in a GeoJSON file containing a rhealpix ID property named 'rhealpix'")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input Rhealpix in GeoJSON"
    )
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of Rhealpix to be uncompacted [0..15]")


    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    
    if resolution < 0 or resolution > 15:
        print(f"Please select a resolution in [0..30] range and try again ")
        return


    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = rhealpixuncompact(rhealpix_dggs,geojson_data,resolution)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"rhealpix_{resolution}_uncompacted.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


import argparse
import pandas as pd
from shapely.geometry import Polygon, box
from shapely.wkt import loads
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.geocode.geocode2geojson import fix_eaggr_wkt
from pyproj import Geod
import fiona
import math
from tqdm import tqdm  # Import tqdm for progress bars
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.geocode.latlon2geocode import latlon2eaggrisea4t
from vgrid.geocode.geocode2geojson import eaggrisea4t2geojson
import geopandas as gpd

import json
import logging
from shapely.geometry import shape, mapping

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize the DGGS model
eaggr_dggs = Eaggr(Model.ISEA4T)

resolution_area = {
    22:1.452590386,
    21:5.810446285,
    20:23.24175429,
    19:92.96709061,
    18:371.8686745,
    17:1487.476609,
    16:5949.921705,
    15:23799.80841,
    14:95200.20638,
    13:380808.5999,
    12:1523297.174,
    11:6091690.181,
    10:24362740.6,
    9:97419500.16,
    8:389424158.4,
    7:1555822168,
    6:6209987762,
    5:24757744008,
    4:98287038022,
    3:3.92934E+11,
    2:1.75574E+12,
    1:6.08082E+12,
    0:2.55613E+13
}

def calculate_bounding_box_area_geod(bbox):
    """
    Calculate the area of a bounding box using WGS84 ellipsoid.
    bbox: tuple of (min_lon, min_lat, max_lon, max_lat)
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    geod = Geod(ellps="WGS84")
    
    # Define the bounding box as a closed loop of points
    lons = [min_lon, max_lon, max_lon, min_lon, min_lon]
    lats = [min_lat, min_lat, max_lat, max_lat, min_lat]
    
    # Calculate the area (negative indicates clockwise orientation)
    area, _ = geod.polygon_area_perimeter(lons, lats)
    
    # Return the absolute area (always positive)
    return abs(area)

def cellid2feature(cell_id):
    eaggr_cell_shape = DggsShape(DggsCell(cell_id), DggsShapeLocation.ONE_FACE)._shape
    cell_to_shp = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell_shape, ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_eaggr_wkt(cell_to_shp)

    # Convert WKT to Shapely geometry
    cell_polygon = loads(cell_to_shp_fixed)
    
    resolution = len(cell_id) - 2
    # Compute centroid
    cell_centroid = cell_polygon.centroid
    center_lat, center_lon = round(cell_centroid.y, 7), round(cell_centroid.x, 7)
    # Compute area using PyProj Geod
    geod = Geod(ellps="WGS84")
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
    # Compute perimeter using PyProj Geod
    edge_len = abs(geod.geometry_area_perimeter(cell_polygon)[1]) / 3  # Perimeter in meters/3

    edge_len_str = f'{round(edge_len, 2)} m'
    cell_area_str = f'{round(cell_area, 2)} m2'

    if cell_area >= 1000_000:
        edge_len_str = f'{round(edge_len / 1000, 2)} km'
        cell_area_str = f'{round(cell_area / (10**6), 2)} km2'

    coordinates = list(cell_polygon.exterior.coords)

    feature = {
        "type": "Feature",
        "geometry": {
            "type": "Polygon",
            "coordinates": [coordinates]  # Wrap the coordinates in an additional list (GeoJSON requirement)
        },
        "properties": {
            "cell_id": cell_id
        }
    }
    
    return feature

def find_appropriate_resolution(bbox_area):
    appropriate_res = 0
    for resolution, avg_area in resolution_area.items():
        if avg_area > bbox_area:
            return resolution


def get_cell_id_covering_bbox_geod(bbox):
    """
    Get the cell ID that covers the bounding box with the smallest resolution.
    bbox: tuple of (min_lon, min_lat, max_lon, max_lat)
    latlon2eaggrisea4t: function to get cell ID
    resolution_area: dictionary mapping resolution to cell area
    """
    # Compute the bounding box area
    bbox_area = calculate_bounding_box_area_geod(bbox)
    # Find the appropriate resolution
    resolution = find_appropriate_resolution(bbox_area) -4 
    if resolution <0:
        resolution = 0
    
    # Calculate the center of the bounding box
    min_lon, min_lat, max_lon, max_lat = bbox
    center_lat = (min_lat + max_lat) / 2.0
    center_lon = (min_lon + max_lon) / 2.0
    
    # Get the cell ID at the resolution
    return latlon2eaggrisea4t(center_lat, center_lon, resolution), resolution


# def get_dggs_cell_siblings(cell_id):

#     cell_siblings = eaggr_dggs.get_dggs_cell_siblings(DggsCell(cell_id))
    
#     return cell_siblings


def get_dggs_cell_children(cell_id, current_resolution, target_resolution):
    """
    Recursively get child cell IDs from current_resolution to target_resolution.
    
    :param cell_id: The initial cell ID at current_resolution.
    :param current_resolution: The current resolution of the given cell ID.
    :param target_resolution: The target resolution to stop the recursion.
    :return: A list of child cell IDs down to the target resolution.
    """
    # Base case: If we've reached the target resolution, return the current cell
    if current_resolution == target_resolution:
        return [cell_id]
    
    # Recursively get the children of the current cell at the next resolution
    children = eaggr_dggs.get_dggs_cell_children(DggsCell(cell_id))
    
    return children


def get_dggs_cell_children_recursive(cell_id, current_resolution, target_resolution):
    """
    Recursively get child cell IDs from current_resolution to target_resolution.
    
    :param cell_id: The initial cell ID at current_resolution.
    :param current_resolution: The current resolution of the given cell ID.
    :param target_resolution: The target resolution to stop the recursion.
    :return: A list of child cell IDs down to the target resolution.
    """
    # Base case: If we've reached the target resolution, return the current cell
    if current_resolution == target_resolution:
        return [cell_id]
    
    # Increment the length of the cell_id (i.e., increase resolution by 1)
    next_resolution = current_resolution + 1
    # Get the child cells at the next resolution
    children = eaggr_dggs.get_dggs_cell_children(DggsCell(cell_id))
    
    # Prepare the list of child cell IDs for the next resolution
    child_cell_ids = []
    for child in children:
        # The child ID at the next resolution is simply the child cell ID
        child_cell_id = child.get_cell_id()
        child_cell_ids.append(child_cell_id)
    
    # Recursively generate child cells for the next resolution
    all_children = []
    for child_cell_id in child_cell_ids:
        all_children.extend(get_dggs_cell_children_recursive(child_cell_id, next_resolution, target_resolution))
    
    return all_children


def main():
    parser = argparse.ArgumentParser(description="Generate DGGS grid within a bounding box and save as GeoJSON.")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Desired resolution (e.g., 0, 1, 2, 3).")
    parser.add_argument(
        "-b", "--bbox", type=float, nargs=4, required=True, 
        help="Bounding box (min_lon, min_lat, max_lon, max_lat)."
    )
    parser.add_argument("-o", "--output", required=True, help="Output GeoJSON path.")
    
    args = parser.parse_args()
    bbox = tuple(args.bbox)
    #  eaggrisea4tgrid -r 14 -b 106.7042196312 10.7741736523 106.7124593773 10.782057247 -o output.geojson
    # Get the bounding box resolution
    resolution = args.resolution
    cell_id, cell_id_res = (get_cell_id_covering_bbox_geod(bbox))
    print(cell_id, cell_id_res)
    siblings = eaggr_dggs.get_dggs_cell_siblings(DggsCell(cell_id))
    siblings.append(DggsCell(cell_id))
    geojson_features = []
    for sibling in siblings:
        # child_cells = get_dggs_cell_children_recursive(cell_id,cell_id_res,resolution)
        child_cells = get_dggs_cell_children_recursive(sibling.get_cell_id(),cell_id_res,cell_id_res + 5)

        # for child_cell in child_cells:
        #     print(child_cell)
        #     geojson_features = []
        
        # Append each child cell as a GeoJSON feature
        for child_cell in child_cells:
            print(f"Processing child cell: {child_cell}")
            # Generate the feature from the child cell ID
            cell_feature = cellid2feature(child_cell)
            geojson_features.append(cell_feature)  # Append the feature to the list
    
    # Save the geojson features to a GeoJSON file
    geojson_output = {
        "type": "FeatureCollection",
        "features": geojson_features  # Append all features here
    }
    
    # Write the GeoJSON to file
    with open(args.output, 'w') as f:
        json.dump(geojson_output, f, indent=2)
    
    print(f"GeoJSON saved to {args.output}")


    # print(cellid2feature(cell_id))
  
if __name__ == "__main__":
    main()

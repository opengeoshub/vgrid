from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from shapely.geometry import Polygon, mapping
from pyproj import Geod
import json, argparse
from vgrid.utils.antimeridian import fix_polygon

isea3h_dggs = Eaggr(Model.ISEA3H)
from pyproj import Geod
geod = Geod(ellps="WGS84")


def cell_to_polygon(isea3h_cell):
    cell_to_shape =  isea3h_dggs.convert_dggs_cell_outline_to_shape_string(isea3h_cell, ShapeStringFormat.WKT)
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


def siblings2geojson(isea3h):
    input_cell = DggsCell(isea3h)
    siblings = isea3h_dggs.get_dggs_cell_siblings(input_cell)
    features = []

    # Process each parent cell
    for sibling in siblings:
        # Convert parent cell to WKT shape string
        cell_polygon = cell_to_polygon(sibling)
        sibling_id = sibling.get_cell_id()
        # if sibling_id.startswith('00') or sibling_id.startswith('03') or sibling_id.startswith('09')\
        # or sibling_id.startswith('14') or sibling_id.startswith('19') or sibling_id.startswith('04') :
        #     cell_polygon = fix_isea3h_antimeridian_cells(cell_polygon)            
   
        cell_centroid = cell_polygon.centroid
        center_lat =  round(cell_centroid.y, 7)
        center_lon = round(cell_centroid.x, 7)
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        avg_edge_len = round(cell_perimeter / 6,2)
        resolution = len(isea3h)-7
        
        # Step 3: Construct the GeoJSON feature
        feature = {
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "isea3h": sibling_id,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_area": cell_area,
                    "avg_edge_len": avg_edge_len,
                    "resolution": resolution
                    }
        }

        features.append(feature)

    # Construct and return the GeoJSON FeatureCollection
    feature_collection = {
        "type": "FeatureCollection",
        "features": features,
    }
    return feature_collection

def siblings2geojson_cli():
    parser = argparse.ArgumentParser(description="Convert OpenEAGGR ISEA3H code to GeoJSON")
    parser.add_argument("isea3h", help="Input ISEA3H code, e.g., isea3h2geojson '07024,0'")
    args = parser.parse_args()
    geojson_data = json.dumps(siblings2geojson(args.isea3h))
    print(geojson_data)
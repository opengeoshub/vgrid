from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.geocode.geocode2geojson import fix_eaggr_wkt
from shapely.wkt import loads
import json
from pyproj import Geod

base_cells = ['00', '01', '02', '03', '04', '06', '06', '07', '08', '09',
              '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
eaggr_dggs = Eaggr(Model.ISEA4T)

def eaggrisea4t2geojson(eaggrisea4t):
    eaggr_cell_shape = DggsShape(DggsCell(eaggrisea4t), DggsShapeLocation.ONE_FACE)._shape
    cell_to_shp = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell_shape, ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_eaggr_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)

    resolution = len(eaggrisea4t) - 2
    # Compute centroid
    cell_centroid = cell_polygon.centroid
    center_lat, center_lon = round(cell_centroid.y, 7), round(cell_centroid.x, 7)
    # Compute area using PyProj Geod
    geod = Geod(ellps="WGS84")
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
    # Compute perimeter using PyProj Geod
    edge_len = abs(geod.geometry_area_perimeter(cell_polygon)[1]) / 3  # Perimeter in meters/ 3

    edge_len_str = f'{round(edge_len, 2)} m'
    cell_area_str = f'{round(cell_area, 2)} m2'

    if cell_area >= 1000_000:
        edge_len_str = f'{round(edge_len / 1000, 2)} km'
        cell_area_str = f'{round(cell_area / (10**6), 2)} km2'

    if cell_polygon:
        coordinates = list(cell_polygon.exterior.coords)
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [coordinates]
            },
            "properties": {
                "eaggr_isea4t": eaggrisea4t,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "cell_area": cell_area_str,
                "edge_len": edge_len_str,
                "resolution": resolution,
            }
        }
        return feature

resolution = 1
# Specify resolution
features = []
feature_collection = None
if (resolution ==0):
    for base_cell in base_cells:
        features.append(eaggrisea4t2geojson(base_cell))
    feature_collection = {
        "type": "FeatureCollection",
        "features": features  # Correctly reference the list of features
    }

    print(json.dumps(feature_collection))
elif resolution == 1:
    for base_cell in base_cells:
        res1_children = children = eaggr_dggs.get_dggs_cell_children(DggsCell((base_cell)))
        for child in res1_children:
            features.append(eaggrisea4t2geojson(child._cell_id))
        feature_collection = {
        "type": "FeatureCollection",
        "features": features  # Correctly reference the list of features
    }
# Save to a GeoJSON file
with open("eaggrisea4t_grid_1.geojson", "w") as f:
    json.dump(feature_collection, f, indent=2)

print("GeoJSON saved as eaggrisea4t_grid_1.geojson")

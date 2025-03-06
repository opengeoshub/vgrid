from shapely.geometry import box, shape, Point, LineString, Polygon, mapping
import argparse
import os, json
from tqdm import tqdm
from vgrid.utils import qtm
from pyproj import Geod
geod = Geod(ellps="WGS84")

p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)


chunk_size = 10_000  # Adjust the chunk size as needed

def chunk_list(data, chunk_size):
    """Yield successive chunks from a list."""
    for i in range(0, len(data), chunk_size):
        yield data[i:i + chunk_size]

# Function to generate grid for Point
def point_to_grid(resolution, point):
    features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    qtm_id = qtm.latlon_to_qtm_id(latitude, longitude, resolution) 
    facet = qtm.qtm_id_to_facet(qtm_id)
    cell_polygon = qtm.constructGeometry(facet)   
    cell_centroid = cell_polygon.centroid
    center_lat =  round(cell_centroid.y, 7)
    center_lon = round(cell_centroid.x, 7)
    cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
    cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
    avg_edge_len = round(cell_perimeter / 3,2)    
    features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                        "qtm": qtm_id,
                        "resolution": resolution,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "avg_edge_len": avg_edge_len,
                        "cell_area": cell_area
                    }
                })
    return {
        "type": "FeatureCollection",
        "features": features,
    }

# Function to generate grid for Polyline
def polyline_to_grid(resolution, geometry):    
    # Extract points from polyline
    if geometry.geom_type == 'LineString':
        # Handle single Polygon as before
        polylines = [geometry]
    elif geometry.geom_type == 'MultiLineString':
        # Handle MultiPolyline: process each polyline separately
        polylines = list(geometry)

    for polyline in polylines:
        # minx, miny, maxx, maxy = polyline.bounds
        # Create a bounding box polygon
        # bbox_poly = box(minx, miny, maxx, maxy)
        levelFacets = {}
        QTMID = {}
        geojson_features = []    
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []

            if lvl == 0:
                initial_facets = [
                    [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                    [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                    [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                    [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                    [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                    [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                    [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                    [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
                ]

                for i, facet in enumerate(initial_facets):
                    QTMID[0].append(str(i + 1))
                    facet_geom = qtm.constructGeometry(facet)
                    cell_centroid = facet_geom.centroid
                    center_lat =  round(cell_centroid.y, 7)
                    center_lon = round(cell_centroid.x, 7)
                    cell_area = round(abs(geod.geometry_area_perimeter(facet_geom)[0]),2)
                    cell_perimeter = abs(geod.geometry_area_perimeter(facet_geom)[1])
                    avg_edge_len = round(cell_perimeter / 3,2)
                    
                    levelFacets[0].append(facet)
                    if shape(facet_geom).intersects(polyline) and resolution == 1 :
                        geojson_features.append({
                            "type": "Feature",
                            "geometry": mapping(facet_geom),
                            "properties": {
                                "qtm": QTMID[0][i],
                                "resolution": resolution,
                                "center_lat": center_lat,
                                "center_lon": center_lon,
                                "avg_edge_len": avg_edge_len,
                                "cell_area": cell_area
                            }
                        })
                        return {
                                "type": "FeatureCollection",
                                "features": geojson_features
                            }              
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if shape(subfacet_geom).intersects(polyline):  # Only keep intersecting facets
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:  # Only store final resolution in GeoJSON
                                cell_centroid = subfacet_geom.centroid
                                center_lat =  round(cell_centroid.y, 7)
                                center_lon = round(cell_centroid.x, 7)
                                cell_area = round(abs(geod.geometry_area_perimeter(subfacet_geom)[0]),2)
                                cell_perimeter = abs(geod.geometry_area_perimeter(subfacet_geom)[1])
                                avg_edge_len = round(cell_perimeter / 3,2)
                                
                                geojson_features.append({
                                    "type": "Feature",
                                    "geometry": mapping(subfacet_geom),
                                    "properties": {
                                        "qtm": new_id,
                                        "resolution": resolution,
                                        "center_lat": center_lat,
                                        "center_lon": center_lon,
                                        "avg_edge_len": avg_edge_len,
                                        "cell_area": cell_area
                                        }
                                })
                            
    return {
            "type": "FeatureCollection",
            "features": geojson_features,
            }


# Function to generate grid for Polygon
def polygon_to_grid(resolution, geometry):
    # Extract points from polyline
    if geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        # Handle MultiPolyline: process each polyline separately
        polygons = list(geometry)

    for polygon in polygons:
        # minx, miny, maxx, maxy = polygon.bounds
        # Create a bounding box polygon
        # bbox_poly = box(minx, miny, maxx, maxy)
        levelFacets = {}
        QTMID = {}
        geojson_features = []    
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []

            if lvl == 0:
                initial_facets = [
                    [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                    [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                    [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                    [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                    [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                    [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                    [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                    [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
                ]

                for i, facet in enumerate(initial_facets):
                    QTMID[0].append(str(i + 1))
                    facet_geom = qtm.constructGeometry(facet)
                    cell_centroid = facet_geom.centroid
                    center_lat =  round(cell_centroid.y, 7)
                    center_lon = round(cell_centroid.x, 7)
                    cell_area = round(abs(geod.geometry_area_perimeter(facet_geom)[0]),2)
                    cell_perimeter = abs(geod.geometry_area_perimeter(facet_geom)[1])
                    avg_edge_len = round(cell_perimeter / 3,2)
                    
                    levelFacets[0].append(facet)
                    if shape(facet_geom).intersects(polygon) and resolution == 1 :
                        geojson_features.append({
                            "type": "Feature",
                            "geometry": mapping(facet_geom),
                            "properties": {
                                "qtm": QTMID[0][i],
                                "resolution": resolution,
                                "center_lat": center_lat,
                                "center_lon": center_lon,
                                "avg_edge_len": avg_edge_len,
                                "cell_area": cell_area
                                }
                        })
                        return {
                                "type": "FeatureCollection",
                                "features": geojson_features
                            }              
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if shape(subfacet_geom).intersects(polygon):  # Only keep intersecting facets
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:  # Only store final resolution in GeoJSON
                                cell_centroid = subfacet_geom.centroid
                                center_lat =  round(cell_centroid.y, 7)
                                center_lon = round(cell_centroid.x, 7)
                                cell_area = round(abs(geod.geometry_area_perimeter(subfacet_geom)[0]),2)
                                cell_perimeter = abs(geod.geometry_area_perimeter(subfacet_geom)[1])
                                avg_edge_len = round(cell_perimeter / 3,2)
                                
                                geojson_features.append({
                                    "type": "Feature",
                                    "geometry": mapping(subfacet_geom),
                                    "properties": {
                                        "qtm": new_id,
                                        "resolution": resolution,
                                        "center_lat": center_lat,
                                        "center_lon": center_lon,
                                        "avg_edge_len": avg_edge_len,
                                        "cell_area": cell_area
                                        }
                                })
                            
    return {
            "type": "FeatureCollection",
            "features": geojson_features,
            }

# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Generate QTM grid for shapes in GeoJSON format")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [1..24]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    
    if resolution < 1 or resolution > 24:
        print(f"Please select a resolution in [1..24] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    # Process GeoJSON features in chunks
    for feature_chunk in tqdm(chunk_list(geojson_data['features'], chunk_size), desc=f"Processing chunks of {chunk_size} features"):
        for feature in feature_chunk:
            if feature['geometry']['type'] in ['Point', 'MultiPoint']:
                coordinates = feature['geometry']['coordinates']
                if feature['geometry']['type'] == 'Point':
                    point = Point(coordinates)
                    point_features = point_to_grid(resolution, point)
                    geojson_features.extend(point_features['features'])

                elif feature['geometry']['type'] == 'MultiPoint':
                    for point_coords in coordinates:
                        point = Point(point_coords)
                        point_features = point_to_grid(resolution, point)
                        geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
                coordinates = feature['geometry']['coordinates']
                if feature['geometry']['type'] == 'LineString':
                    polyline = LineString(coordinates)
                    polyline_features = polyline_to_grid(resolution, polyline)
                    geojson_features.extend(polyline_features['features'])

                elif feature['geometry']['type'] == 'MultiLineString':
                    for line_coords in coordinates:
                        polyline = LineString(line_coords)
                        polyline_features = polyline_to_grid(resolution, polyline)
                        geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
                coordinates = feature['geometry']['coordinates']

                if feature['geometry']['type'] == 'Polygon':
                    exterior_ring = coordinates[0]
                    interior_rings = coordinates[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = polygon_to_grid(resolution, polygon)
                    geojson_features.extend(polygon_features['features'])

                elif feature['geometry']['type'] == 'MultiPolygon':
                    for sub_polygon_coords in coordinates:
                        exterior_ring = sub_polygon_coords[0]
                        interior_rings = sub_polygon_coords[1:]
                        polygon = Polygon(exterior_ring, interior_rings)
                        polygon_features = polygon_to_grid(resolution, polygon)
                        geojson_features.extend(polygon_features['features'])

   
    # Save the results to GeoJSON
    geojson_path = f"geojson2qtm_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()

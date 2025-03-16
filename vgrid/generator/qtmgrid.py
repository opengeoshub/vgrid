from shapely.geometry import shape, Polygon, mapping
import argparse
import json
from vgrid.utils import qtm
from pyproj import Geod
geod = Geod(ellps="WGS84")
from vgrid.generator.settings import max_cells

p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)


def generate_grid(resolution):
    """Generates a Dutton QTM grid at a specific resolution and saves it as GeoJSON."""
    
    levelFacets = {}
    QTMID = {}

    for lvl in range(resolution):
        levelFacets[lvl] = []
        QTMID[lvl] = []
        geojson_features = []  # Store GeoJSON features separately

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
                facet_geom = qtm.constructGeometry(facet)
                QTMID[0].append(str(i + 1))
                levelFacets[0].append(facet)
                
                cell_centroid = facet_geom.centroid
                center_lat =  round(cell_centroid.y, 7)
                center_lon = round(cell_centroid.x, 7)
                cell_area = round(abs(geod.geometry_area_perimeter(facet_geom)[0]),2)
                cell_perimeter = abs(geod.geometry_area_perimeter(facet_geom)[1])
                avg_edge_len = round(cell_perimeter / 3,2)                
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
                
        else:
            for i, pf in enumerate(levelFacets[lvl - 1]):
                subdivided_facets = qtm.divideFacet(pf)
                for j, subfacet in enumerate(subdivided_facets):
                    new_id = QTMID[lvl - 1][i] + str(j)
                    QTMID[lvl].append(new_id)                    
                    levelFacets[lvl].append(subfacet)
                    if lvl == resolution -1:
                        subfacet_geom= qtm.constructGeometry(subfacet)
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
        "features": geojson_features
    }

def generate_grid_within_bbox(resolution, bbox):
    """Generates a Dutton QTM grid at a specific resolution within a bounding box and saves it as GeoJSON."""
    levelFacets = {}
    QTMID = {}
    geojson_features = []

    # Convert bbox to Polygon
    bbox_poly = Polygon([
        (bbox[0], bbox[1]),  # min_lon, min_lat
        (bbox[2], bbox[1]),  # max_lon, min_lat
        (bbox[2], bbox[3]),  # max_lon, max_lat
        (bbox[0], bbox[3]),  # min_lon, max_lat
        (bbox[0], bbox[1])   # Close the polygon
    ])

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
                if shape(facet_geom).intersects(bbox_poly) and resolution == 1 :
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
                    if shape(subfacet_geom).intersects(bbox_poly):  # Only keep intersecting facets
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
        "features": geojson_features
    }
    
def main():
    parser = argparse.ArgumentParser(description='Generate QTM grid.')
    parser.add_argument('-r', '--resolution', required=True, type=int, help='Resolution [1..24] to generate.')
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    
    if resolution < 0 or resolution > 24:
        print(f"Please select a resolution in [0..15] range and try again ")
        return

    if bbox == [-180, -90, 180, 90]:
        geojson_features = generate_grid(resolution)
        # Define the GeoJSON file path
        geojson_path = f"qtm_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    
    else:
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"qtm_grid_{resolution}_bbox.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")

if __name__ == '__main__':
    main()

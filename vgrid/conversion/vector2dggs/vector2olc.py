import os
import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon, shape
import pandas as pd
import geopandas as gpd
from vgrid.utils import olc
from vgrid.generator.olcgrid import generate_grid, refine_cell
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.dggscompact import olccompact
from vgrid.utils.geometry import check_predicate

def validate_olc_resolution(resolution):
    """
    Validate that OLC resolution is in the valid range [2,4,6,8,10,11,12,13,14,15].
    
    Args:
        resolution: Resolution value to validate
        
    Returns:
        int: Validated resolution value
        
    Raises:
        ValueError: If resolution is not in range [2,4,6,8,10,11,12,13,14,15]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(f"Resolution must be an integer, got {type(resolution).__name__}")
    
    if resolution not in [2,4,6,8,10,11,12,13,14,15]:
        raise ValueError(f"Resolution must be in [2,4,6,8,10,11,12,13,14,15], got {resolution}")
    
    return resolution


def point2olc(resolution, point, feature_properties=None, predicate = None, compact= False, topology=False, include_properties=True):
    olc_features = []
    olc_id = olc.encode(point.y, point.x, resolution)
    coord = olc.decode(olc_id)
    if coord:
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi
        cell_polygon = Polygon([
            [min_lon, min_lat],
            [max_lon, min_lat],
            [max_lon, max_lat],
            [min_lon, max_lat],
            [min_lon, min_lat]
        ])
        olc_feature = graticule_dggs_to_feature("olc", olc_id, resolution, cell_polygon)
        if include_properties and feature_properties:
            olc_feature["properties"].update(feature_properties)
        olc_features.append(olc_feature)
    return olc_features

def polyline2olc(resolution, feature, feature_properties=None, predicate = None, compact= False, topology=False, include_properties=True):
    olc_features = []
    if feature.geom_type in ('LineString'):
        polylines = [feature]
    elif feature.geom_type in ('MultiLineString'):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        base_resolution = 2
        base_cells = generate_grid(base_resolution, verbose = False)
        seed_cells = []
        for base_cell in base_cells["features"]:
            base_cell_poly = Polygon(base_cell["geometry"]["coordinates"][0])
            if polyline.intersects(base_cell_poly):
                seed_cells.append(base_cell)
        refined_features = []
        for seed_cell in seed_cells:
            seed_cell_poly = Polygon(seed_cell["geometry"]["coordinates"][0])
            if seed_cell_poly.contains(polyline) and resolution == base_resolution:
                refined_features.append(seed_cell)
            else:
                refined_features.extend(
                    refine_cell(seed_cell_poly.bounds, base_resolution, resolution, polyline)
                )
        resolution_features = [
            refined_feature for refined_feature in refined_features if refined_feature["properties"]["resolution"] == resolution
        ]
        seen_olc_codes = set()
        for resolution_feature in resolution_features:
            olc_id = resolution_feature["properties"]["olc"]
            if olc_id not in seen_olc_codes:
                if include_properties and feature_properties:
                    resolution_feature["properties"].update(feature_properties)
                olc_features.append(resolution_feature)
                seen_olc_codes.add(olc_id)
    olc_geojson = {
        "type": "FeatureCollection",
        "features": olc_features
    }
    if compact:
        return olccompact(olc_geojson)["features"]
    return olc_features


def polygon2olc(resolution, feature, feature_properties=None, predicate = None, compact= False, topology=False, include_properties=True):
    olc_features = []
    if feature.geom_type in ('Polygon'):
        polygons = [feature]
    elif feature.geom_type in ('MultiPolygon'):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:  
        base_resolution = 2
        base_cells = generate_grid(base_resolution, verbose = False )
        seed_cells = []
        for base_cell in base_cells["features"]:
            base_cell_poly = Polygon(base_cell["geometry"]["coordinates"][0])
            if polygon.intersects(base_cell_poly):
                seed_cells.append(base_cell)
        refined_features = []
        for seed_cell in seed_cells:
            seed_cell_poly = Polygon(seed_cell["geometry"]["coordinates"][0])
            if seed_cell_poly.contains(polygon) and resolution == base_resolution:
                refined_features.append(seed_cell)
            else:
                refined_features.extend(
                    refine_cell(seed_cell_poly.bounds, base_resolution, resolution, polygon)    
                )
        resolution_features = [
            refined_feature for refined_feature in refined_features if refined_feature["properties"]["resolution"] == resolution
        ]
        seen_olc_codes = set()
        for resolution_feature in resolution_features:
            olc_id = resolution_feature["properties"]["olc"]
            if olc_id not in seen_olc_codes:
                cell_geom = shape(resolution_feature["geometry"])
                if not check_predicate(cell_geom, polygon, predicate):
                    continue
                if include_properties and feature_properties:
                    resolution_feature["properties"].update(feature_properties)
                olc_features.append(resolution_feature)
                seen_olc_codes.add(olc_id)
    olc_geojson = {
        "type": "FeatureCollection",
        "features": olc_features
    }
    if compact:
        return olccompact(olc_geojson)["features"]
    return olc_features

def geometry2olc(geometries, resolution, properties_list=None, predicate=None, compact=False, topology=False, include_properties=True):
    resolution = validate_olc_resolution(resolution)
    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]
    
    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]


    olc_features = []
    for idx, geom in tqdm(enumerate(geometries), desc="Processing features"):
        props = properties_list[idx] if properties_list and idx < len(properties_list) else {}
        if not include_properties:
            props = {}
        if geom.geom_type == 'Point':
            olc_features.extend(point2olc(resolution, geom, props, predicate, compact, topology, include_properties))
        elif geom.geom_type == 'MultiPoint':
            for pt in geom.geoms:
                olc_features.extend(point2olc(resolution, pt, props, predicate, compact, topology, include_properties))
        elif geom.geom_type in ('LineString', 'MultiLineString'):
            olc_features.extend(polyline2olc(resolution, geom, props, predicate, compact, topology, include_properties))
        elif geom.geom_type in ('Polygon', 'MultiPolygon'):
            olc_features.extend(polygon2olc(resolution, geom, props, predicate, compact, topology, include_properties))
    return {"type": "FeatureCollection", "features": olc_features}

def dataframe2olc(df, resolution, predicate=None, compact=False, topology=False, include_properties=True):
    geometries = []
    properties_list = []
    for idx, row in df.iterrows():
        geom = row.geometry if 'geometry' in row else row['geometry']
        if geom is not None:
            geometries.append(geom)
            props = row.to_dict()
            if 'geometry' in props:
                del props['geometry']
            properties_list.append(props)
    return geometry2olc(geometries, resolution, properties_list, predicate, compact, topology, include_properties)

def geodataframe2olc(gdf, resolution, predicate=None, compact=False, topology=False, include_properties=True):
    geometries = []
    properties_list = []
    for idx, row in gdf.iterrows():
        geom = row.geometry
        if geom is not None:
            geometries.append(geom)
            props = row.to_dict()
            if 'geometry' in props:
                del props['geometry']
            properties_list.append(props)
    return geometry2olc(geometries, resolution, properties_list, predicate, compact, topology, include_properties)

def vector2olc(data, resolution, predicate=None, compact=False, topology=False, output_format='geojson', output_path=None, include_properties=True, **kwargs):
    if hasattr(data, 'geometry') and hasattr(data, 'columns'):
        result = geodataframe2olc(data, resolution, predicate, compact, topology, include_properties)
    elif isinstance(data, pd.DataFrame):
        result = dataframe2olc(data, resolution, predicate, compact, topology, include_properties)
    elif hasattr(data, 'geom_type') or (isinstance(data, list) and len(data) > 0 and hasattr(data[0], 'geom_type')):
        result = geometry2olc(data, resolution, None, predicate, compact, topology, include_properties)
    elif isinstance(data, dict) and 'type' in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data['features'])
            result = geodataframe2olc(gdf, resolution, predicate, compact, topology, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2olc(gdf, resolution, predicate, compact, topology, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")
    return convert_to_output_format(result, output_format, output_path)

def convert_to_output_format(result, output_format, output_path=None):
    gdf = gpd.GeoDataFrame.from_features(result['features'])
    gdf.set_crs(epsg=4326, inplace=True)
    if output_format.lower() == 'geojson':
        if output_path:
            with open(output_path, 'w') as f:
                json.dump(result, f, indent=2)
            return output_path
        else:
            return result
    elif output_format.lower() == 'gpkg':
        if output_path:
            gdf.to_file(output_path, driver='GPKG')
            return output_path
        else:
            gdf.to_file("vector2olc.gpkg", driver='GPKG')
            return "vector2olc.gpkg"
    elif output_format.lower() == 'parquet':
        if output_path:
            gdf.to_parquet(output_path, index=False)
            return output_path
        else:
            return gdf.to_parquet(index=False)
    elif output_format.lower() == 'csv':
        if output_path:
            gdf.to_csv(output_path, index=False)
            return output_path
        else:
            return gdf.to_csv(index=False)
    elif output_format.lower() == 'shapefile':
        if output_path:
            gdf.to_file(output_path, driver="ESRI Shapefile")
            return output_path
        else:
            gdf.to_file("vector2olc.shp", driver="ESRI Shapefile")
            return "vector2olc.shp"
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile")

def vector2olc_cli():
    parser = argparse.ArgumentParser(description='Convert vector data to OLC grid cells')
    parser.add_argument('-i', '--input', help='Input file path, URL')
    parser.add_argument('-r', '--resolution', type=int, choices=[2, 4, 6, 8, 10, 11, 12, 13, 14, 15], metavar='[2,4,6,8,10,11,12,13,14,15]', help='OLC resolution (see OLC spec)')
    parser.add_argument('-c', '--compact', action='store_true', help='Enable OLC compact mode for polygons')
    parser.add_argument('-p', '--predicate', 
                       choices=['intersect', 'within', 'centroid_within', 'largest_overlap'],
                       help='Spatial predicate: intersect, within, centroid_within, largest_overlap for polygons')
    parser.add_argument('-t', '--topology', action='store_true', help='Enable topology preserving mode')
    parser.add_argument('-np', '-no-props', dest='include_properties', action='store_false', help="Do not include original feature properties.")
    parser.add_argument('-f', '--format', default='geojson', choices=['geojson', 'gpkg', 'parquet', 'csv', 'shapefile'], help='Output format (default: geojson)')
    parser.add_argument('-o', '--output', help='Output file path (optional)')
    args = parser.parse_args()
    if args.resolution is not None:
        try:
            args.resolution = validate_olc_resolution(args.resolution)
        except (ValueError, TypeError) as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
    data = args.input
    output_path = args.output
    if not output_path and args.format in ['geojson', 'gpkg', 'parquet', 'csv', 'shapefile']:
        extensions = {
            'geojson': '.geojson',
            'gpkg': '.gpkg',
            'parquet': '.parquet',
            'csv': '.csv',
            'shapefile': '.shp'
        }
        output_path = f"vector2olc{extensions.get(args.format, '')}"
    try:
        result = vector2olc(
            data,
            args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            output_format=args.format,
            output_path=output_path,
            include_properties=args.include_properties
        )
        if output_path:
            print(f"Output saved to {output_path}")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    vector2olc_cli() 
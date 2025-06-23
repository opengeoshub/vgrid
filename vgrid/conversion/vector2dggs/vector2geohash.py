import os
import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import geohash
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.generator.geohashgrid import initial_geohashes, geohash_to_polygon, expand_geohash_bbox
from vgrid.conversion.dggscompact import geohashcompact

# --- Geohash Conversion Helpers (adapted from geojson2geohash.py) ---
def point_to_geohash(resolution, point, feature_properties=None):
    geohash_features = []
    longitude = point.x
    latitude = point.y
    geohash_id = geohash.encode(latitude, longitude, resolution)
    bbox = geohash.bbox(geohash_id)
    if bbox:
        min_lat, min_lon = bbox['s'], bbox['w']
        max_lat, max_lon = bbox['n'], bbox['e']
        cell_polygon = Polygon([
            [min_lon, min_lat],
            [max_lon, min_lat],
            [max_lon, max_lat],
            [min_lon, max_lat],
            [min_lon, min_lat]
        ])
        geohash_feature = graticule_dggs_to_feature("geohash", geohash_id, resolution, cell_polygon)
        if feature_properties:
            geohash_feature["properties"].update(feature_properties)
        geohash_features.append(geohash_feature)
    return geohash_features

def poly_to_geohash(resolution, geometry, feature_properties=None, compact=False):
    geohash_features = []
    if geometry.geom_type in ('LineString', 'Polygon'):
        polys = [geometry]
    elif geometry.geom_type in ('MultiLineString', 'MultiPolygon'):
        polys = list(geometry)
    else:
        return []
    for poly in polys:
        intersected_geohashes = {gh for gh in initial_geohashes if geohash_to_polygon(gh).intersects(poly)}
        geohashes_bbox = set()
        for gh in intersected_geohashes:
            expand_geohash_bbox(gh, resolution, geohashes_bbox, poly)
        for gh in geohashes_bbox:
            cell_polygon = geohash_to_polygon(gh)
            geohash_feature = graticule_dggs_to_feature("geohash", gh, resolution, cell_polygon)
            if feature_properties:
                geohash_feature["properties"].update(feature_properties)
            geohash_features.append(geohash_feature)
    geohash_geojson = {
        "type": "FeatureCollection",
        "features": geohash_features
    }
    if compact:
        return geohashcompact(geohash_geojson)["features"]
    return geohash_features

# --- Main geometry conversion ---
def geometry2geohash(geometries, resolution, properties_list=None, compact=False, include_properties=True):
    geohash_features = []
    for idx, geom in enumerate(geometries):
        props = properties_list[idx] if properties_list and idx < len(properties_list) else {}
        if not include_properties:
            props = {}
        if geom.geom_type == 'Point':
            geohash_features.extend(point_to_geohash(resolution, geom, props))
        elif geom.geom_type == 'MultiPoint':
            for pt in geom.geoms:
                geohash_features.extend(point_to_geohash(resolution, pt, props))
        elif geom.geom_type in ('LineString', 'MultiLineString'):
            geohash_features.extend(poly_to_geohash(resolution, geom, props, compact))
        elif geom.geom_type in ('Polygon', 'MultiPolygon'):
            geohash_features.extend(poly_to_geohash(resolution, geom, props, compact))
    return {"type": "FeatureCollection", "features": geohash_features}

# --- DataFrame/GeoDataFrame conversion ---
def dataframe2geohash(df, resolution, compact=False, include_properties=True):
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
    return geometry2geohash(geometries, resolution, properties_list, compact, include_properties)

def geodataframe2geohash(gdf, resolution, compact=False, include_properties=True):
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
    return geometry2geohash(geometries, resolution, properties_list, compact, include_properties)

# --- Main vector2geohash function ---
def vector2geohash(data, resolution, compact=False, output_format='geojson', output_path=None, include_properties=True, **kwargs):
    """
    Convert vector data to Geohash grid cells from various input formats.
    Args:
        data: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
        resolution (int): Geohash resolution [1..10]
        compact (bool): Enable Geohash compact mode for polygons (default: False)
        output_format (str): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str): Output file path (optional)
        include_properties (bool): If False, do not include original feature properties. (default: True)
        **kwargs: Additional arguments passed to geopandas read functions
    Returns:
        dict or str: Output in the specified format
    """
    if not isinstance(resolution, int) or resolution < 1 or resolution > 10:
        raise ValueError(f"Resolution must be in range [1..10], got {resolution}")
    # Process input data directly
    if hasattr(data, 'geometry') and hasattr(data, 'columns'):
        # GeoDataFrame
        result = geodataframe2geohash(data, resolution, compact, include_properties)
    elif isinstance(data, pd.DataFrame):
        # DataFrame with geometry column
        result = dataframe2geohash(data, resolution, compact, include_properties)
    elif hasattr(data, 'geom_type') or (isinstance(data, list) and len(data) > 0 and hasattr(data[0], 'geom_type')):
        # Shapely geometry objects
        result = geometry2geohash(data, resolution, None, compact, include_properties)
    elif isinstance(data, dict) and 'type' in data:
        # GeoJSON dict
        try:
            gdf = gpd.GeoDataFrame.from_features(data['features'])
            result = geodataframe2geohash(gdf, resolution, compact, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        # File path or URL
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2geohash(gdf, resolution, compact, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")
    return convert_to_output_format(result, output_format, output_path)

# --- Output format conversion ---
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
            gdf.to_file("vector2geohash.gpkg", driver='GPKG')
            return "vector2geohash.gpkg"
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
            gdf.to_file("vector2geohash.shp", driver="ESRI Shapefile")
            return "vector2geohash.shp"
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile")

# --- CLI ---
def vector2geohash_cli():
    parser = argparse.ArgumentParser(description='Convert vector data to Geohash grid cells')
    parser.add_argument('-i', '--input', help='Input file path, URL')
    parser.add_argument('-r', '--resolution', type=int, choices=range(1, 11), metavar='[1-10]', help='Geohash resolution [1..10] (1=coarsest, 10=finest)')
    parser.add_argument('-c', '--compact', action='store_true', help='Enable Geohash compact mode for polygons')
    parser.add_argument('-np', '-no-props', dest='include_properties', action='store_false', help="Do not include original feature properties.")
    parser.add_argument('-f', '--format', default='geojson', choices=['geojson', 'gpkg', 'parquet', 'csv', 'shapefile'], help='Output format (default: geojson)')
    parser.add_argument('-o', '--output', help='Output file path (optional)')
    args = parser.parse_args()
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
        output_path = f"vector2geohash{extensions.get(args.format, '')}"
    try:
        result = vector2geohash(
            data,
            args.resolution,
            compact=args.compact,
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
    vector2geohash_cli() 
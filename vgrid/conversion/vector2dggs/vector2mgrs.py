import os
import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import mgrs
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.latlon2dggs import latlon2mgrs
from vgrid.conversion.dggs2geojson import mgrs2geojson

def point_to_mgrs(resolution, point, feature_properties=None):
    mgrs_features = []
    latitude, longitude = point.y, point.x
    mgrs_id = latlon2mgrs(latitude, longitude, resolution)
    mgrs_feature = mgrs2geojson(mgrs_id)
    if mgrs_feature:
        if feature_properties:
            mgrs_feature["properties"].update(feature_properties)
        mgrs_features.append(mgrs_feature)
    return mgrs_features

def poly_to_mgrs(resolution, geometry, feature_properties=None):
    # Placeholder: MGRS polyline/polygon support not implemented in geojson2mgrs.py
    return []

def geometry2mgrs(geometries, resolution, properties_list=None, compact=False, include_properties=True):
    mgrs_features = []
    for idx, geom in enumerate(geometries):
        props = properties_list[idx] if properties_list and idx < len(properties_list) else {}
        if not include_properties:
            props = {}
        if geom.geom_type == 'Point':
            mgrs_features.extend(point_to_mgrs(resolution, geom, props))
        elif geom.geom_type == 'MultiPoint':
            for pt in geom.geoms:
                mgrs_features.extend(point_to_mgrs(resolution, pt, props))
        elif geom.geom_type in ('LineString', 'MultiLineString', 'Polygon', 'MultiPolygon'):
            mgrs_features.extend(poly_to_mgrs(resolution, geom, props))
    return {"type": "FeatureCollection", "features": mgrs_features}

def dataframe2mgrs(df, resolution, compact=False, include_properties=True):
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
    return geometry2mgrs(geometries, resolution, properties_list, compact, include_properties)

def geodataframe2mgrs(gdf, resolution, compact=False, include_properties=True):
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
    return geometry2mgrs(geometries, resolution, properties_list, compact, include_properties)

def vector2mgrs(data, resolution, compact=False, output_format='geojson', output_path=None, include_properties=True, **kwargs):
    if not isinstance(resolution, int) or resolution < 0 or resolution > 5:
        raise ValueError(f"Resolution must be in range [0..5], got {resolution}")
    if hasattr(data, 'geometry') and hasattr(data, 'columns'):
        result = geodataframe2mgrs(data, resolution, compact, include_properties)
    elif isinstance(data, pd.DataFrame):
        result = dataframe2mgrs(data, resolution, compact, include_properties)
    elif hasattr(data, 'geom_type') or (isinstance(data, list) and len(data) > 0 and hasattr(data[0], 'geom_type')):
        result = geometry2mgrs(data, resolution, None, compact, include_properties)
    elif isinstance(data, dict) and 'type' in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data['features'])
            result = geodataframe2mgrs(gdf, resolution, compact, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2mgrs(gdf, resolution, compact, include_properties)
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
            gdf.to_file("vector2mgrs.gpkg", driver='GPKG')
            return "vector2mgrs.gpkg"
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
            gdf.to_file("vector2mgrs.shp", driver="ESRI Shapefile")
            return "vector2mgrs.shp"
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile")

def vector2mgrs_cli():
    parser = argparse.ArgumentParser(description='Convert vector data to MGRS grid cells')
    parser.add_argument('-i', '--input', help='Input file path, URL')
    parser.add_argument('-r', '--resolution', type=int, choices=range(0, 6), metavar='[0-5]', help='MGRS resolution [0..5] (0=coarsest, 5=finest)')
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
        output_path = f"vector2mgrs{extensions.get(args.format, '')}"
    try:
        result = vector2mgrs(
            data,
            args.resolution,
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
    vector2mgrs_cli() 
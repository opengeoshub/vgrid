import os
import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import s2
from vgrid.generator.s2grid import s2_cell_to_polygon
from vgrid.generator.settings import geodesic_dggs_to_feature

# --- S2 Conversion Helpers (adapted from geojson2s2.py) ---
def point_to_s2(resolution, point, feature_properties=None):
    s2_features = []
    latitude = point.y
    longitude = point.x
    lat_lng = s2.LatLng.from_degrees(latitude, longitude)
    cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
    cell_id = cell_id_max_res.parent(resolution)
    s2_cell = s2.Cell(cell_id)
    cell_token = s2.CellId.to_token(s2_cell.id())
    if s2_cell:
        cell_polygon = s2_cell_to_polygon(cell_id)
        resolution = cell_id.level()
        num_edges = 4
        s2_feature = geodesic_dggs_to_feature("s2", cell_token, resolution, cell_polygon, num_edges)
        if feature_properties:
            s2_feature["properties"].update(feature_properties)
        s2_features.append(s2_feature)
    return s2_features

def poly_to_s2(resolution, geometry, feature_properties=None, compact=None):
    s2_features = []
    if geometry.geom_type in ('LineString', 'Polygon'):
        polys = [geometry]
    elif geometry.geom_type in ('MultiLineString', 'MultiPolygon'):
        polys = list(geometry)
    else:
        return []
    for poly in polys:
        min_lng, min_lat, max_lng, max_lat = poly.bounds
        level = resolution
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng)
        )
        covering = coverer.get_covering(region)
        cell_ids = covering
        if compact:
            covering = s2.CellUnion(covering)
            covering.normalize()
            cell_ids = covering.cell_ids()
        for cell_id in cell_ids:
            cell_polygon = s2_cell_to_polygon(cell_id)
            if cell_polygon.intersects(poly):
                cell_token = s2.CellId.to_token(cell_id)
                cell_resolution = cell_id.level()
                num_edges = 4
                s2_feature = geodesic_dggs_to_feature("s2", cell_token, cell_resolution, cell_polygon, num_edges)
                if feature_properties:
                    s2_feature["properties"].update(feature_properties)
                s2_features.append(s2_feature)
    return s2_features

# --- Main geometry conversion ---
def geometry2s2(geometries, resolution, properties_list=None, compact=False, include_properties=True):
    s2_features = []
    for idx, geom in enumerate(geometries):
        props = properties_list[idx] if properties_list and idx < len(properties_list) else {}
        if not include_properties:
            props = {}
        if geom.geom_type == 'Point':
            s2_features.extend(point_to_s2(resolution, geom, props))
        elif geom.geom_type == 'MultiPoint':
            for pt in geom.geoms:
                s2_features.extend(point_to_s2(resolution, pt, props))
        elif geom.geom_type in ('LineString', 'MultiLineString'):
            s2_features.extend(poly_to_s2(resolution, geom, props, compact))
        elif geom.geom_type in ('Polygon', 'MultiPolygon'):
            s2_features.extend(poly_to_s2(resolution, geom, props, compact))
    return {"type": "FeatureCollection", "features": s2_features}

# --- DataFrame/GeoDataFrame conversion ---
def dataframe2s2(df, resolution, compact=False, include_properties=True):
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
    return geometry2s2(geometries, resolution, properties_list, compact, include_properties)

def geodataframe2s2(gdf, resolution, compact=False, include_properties=True):
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
    return geometry2s2(geometries, resolution, properties_list, compact, include_properties)

# --- Main vector2s2 function ---
def vector2s2(data, resolution, compact=False, output_format='geojson', output_path=None, include_properties=True, **kwargs):
    """
    Convert vector data to S2 grid cells from various input formats.
    Args:
        data: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
        resolution (int): S2 resolution [0..30]
        compact (bool): Enable S2 compact mode for polygons (default: False)
        output_format (str): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str): Output file path (optional)
        include_properties (bool): If False, do not include original feature properties. (default: True)
        **kwargs: Additional arguments passed to geopandas read functions
    Returns:
        dict or str: Output in the specified format
    """
    if not isinstance(resolution, int) or resolution < 0 or resolution > 30:
        raise ValueError(f"Resolution must be in range [0..30], got {resolution}")
    # Process input data directly
    if hasattr(data, 'geometry') and hasattr(data, 'columns'):
        # GeoDataFrame
        result = geodataframe2s2(data, resolution, compact, include_properties)
    elif isinstance(data, pd.DataFrame):
        # DataFrame with geometry column
        result = dataframe2s2(data, resolution, compact, include_properties)
    elif hasattr(data, 'geom_type') or (isinstance(data, list) and len(data) > 0 and hasattr(data[0], 'geom_type')):
        # Shapely geometry objects
        result = geometry2s2(data, resolution, None, compact, include_properties)
    elif isinstance(data, dict) and 'type' in data:
        # GeoJSON dict
        try:
            gdf = gpd.GeoDataFrame.from_features(data['features'])
            result = geodataframe2s2(gdf, resolution, compact, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        # File path or URL
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2s2(gdf, resolution, compact, include_properties)
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
            gdf.to_file("vector2s2.gpkg", driver='GPKG')
            return "vector2s2.gpkg"
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
            gdf.to_file("vector2s2.shp", driver="ESRI Shapefile")
            return "vector2s2.shp"
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile")

# --- CLI ---
def vector2s2_cli():
    parser = argparse.ArgumentParser(description='Convert vector data to S2 grid cells')
    parser.add_argument('-i', '--input', help='Input file path, URL')
    parser.add_argument('-r', '--resolution', type=int, choices=range(31), metavar='[0-30]', help='S2 resolution [0..30] (0=coarsest, 30=finest)')
    parser.add_argument('-c', '--compact', action='store_true', help='Enable S2 compact mode for polygons')
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
        output_path = f"vector2s2{extensions.get(args.format, '')}"
    try:
        result = vector2s2(
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
    vector2s2_cli() 
import os
import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon, box
import pandas as pd
import geopandas as gpd
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from vgrid.utils.rhealpixdggs.conversion import compress_order_cells, get_finest_containing_cell
from vgrid.generator.rhealpixgrid import fix_rhealpix_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.conversion.dggscompact import rhealpix_compact

def rhealpix_cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = fix_rhealpix_antimeridian_cells(vertices)
    return Polygon(vertices)

def point_to_rhealpix(rhealpix_dggs, resolution, point, feature_properties=None):
    rhealpix_features = []
    seed_cell = rhealpix_dggs.cell_from_point(resolution, (point.x, point.y), plane=False)
    seed_cell_polygon = rhealpix_cell_to_polygon(seed_cell)
    if seed_cell_polygon:
        seed_cell_id = str(seed_cell)
        num_edges = 4
        if seed_cell.ellipsoidal_shape() == 'dart':
            num_edges = 3
        rhealpix_feature = geodesic_dggs_to_feature("rhealpix", seed_cell_id, resolution, seed_cell_polygon, num_edges)
        if feature_properties:
            rhealpix_feature["properties"].update(feature_properties)
        rhealpix_features.append(rhealpix_feature)
    return rhealpix_features

def poly_to_rhealpix(rhealpix_dggs, resolution, geometry, feature_properties=None, compact=False):
    rhealpix_features = []
    if geometry.geom_type in ('LineString', 'Polygon'):
        polys = [geometry]
    elif geometry.geom_type in ('MultiLineString', 'MultiPolygon'):
        polys = list(geometry)
    else:
        return []
    for poly in polys:
        minx, miny, maxx, maxy = poly.bounds
        bbox_polygon = box(minx, miny, maxx, maxy)
        bbox_center_lon = bbox_polygon.centroid.x
        bbox_center_lat = bbox_polygon.centroid.y
        seed_point = (bbox_center_lon, bbox_center_lat)
        seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
        seed_cell_id = str(seed_cell)
        seed_cell_polygon = rhealpix_cell_to_polygon(seed_cell)
        if seed_cell_polygon.contains(bbox_polygon):
            num_edges = 4
            if seed_cell.ellipsoidal_shape() == 'dart':
                num_edges = 3
            cell_resolution = resolution
            rhealpix_feature = geodesic_dggs_to_feature("rhealpix", seed_cell_id, cell_resolution, seed_cell_polygon, num_edges)
            if feature_properties:
                rhealpix_feature["properties"].update(feature_properties)
            rhealpix_features.append(rhealpix_feature)
            return rhealpix_features
        else:
            covered_cells = set()
            queue = [seed_cell]
            while queue:
                current_cell = queue.pop()
                current_cell_id = str(current_cell)
                if current_cell_id in covered_cells:
                    continue
                covered_cells.add(current_cell_id)
                cell_polygon = rhealpix_cell_to_polygon(current_cell)
                if not cell_polygon.intersects(bbox_polygon):
                    continue
                neighbors = current_cell.neighbors(plane=False)
                for _, neighbor in neighbors.items():
                    neighbor_id = str(neighbor)
                    if neighbor_id not in covered_cells:
                        queue.append(neighbor)
            if compact:
                covered_cells = rhealpix_compact(rhealpix_dggs, covered_cells)
            for cell_id in covered_cells:
                rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
                rhelpix_cell = rhealpix_dggs.cell(rhealpix_uids)
                cell_resolution = rhelpix_cell.resolution
                cell_polygon = rhealpix_cell_to_polygon(rhelpix_cell)
                if cell_polygon.intersects(poly):
                    num_edges = 4
                    if seed_cell.ellipsoidal_shape() == 'dart':
                        num_edges = 3
                    rhealpix_feature = geodesic_dggs_to_feature("rhealpix", str(cell_id), cell_resolution, cell_polygon, num_edges)
                    if feature_properties:
                        rhealpix_feature["properties"].update(feature_properties)
                    rhealpix_features.append(rhealpix_feature)
    return rhealpix_features

def geometry2rhealpix(geometries, resolution, properties_list=None, compact=False, include_properties=True):
    rhealpix_dggs = RHEALPixDGGS()
    rhealpix_features = []
    for idx, geom in enumerate(geometries):
        props = properties_list[idx] if properties_list and idx < len(properties_list) else {}
        if not include_properties:
            props = {}
        if geom.geom_type == 'Point':
            rhealpix_features.extend(point_to_rhealpix(rhealpix_dggs, resolution, geom, props))
        elif geom.geom_type == 'MultiPoint':
            for pt in geom.geoms:
                rhealpix_features.extend(point_to_rhealpix(rhealpix_dggs, resolution, pt, props))
        elif geom.geom_type in ('LineString', 'MultiLineString'):
            rhealpix_features.extend(poly_to_rhealpix(rhealpix_dggs, resolution, geom, props, compact))
        elif geom.geom_type in ('Polygon', 'MultiPolygon'):
            rhealpix_features.extend(poly_to_rhealpix(rhealpix_dggs, resolution, geom, props, compact))
    return {"type": "FeatureCollection", "features": rhealpix_features}

def dataframe2rhealpix(df, resolution, compact=False, include_properties=True):
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
    return geometry2rhealpix(geometries, resolution, properties_list, compact, include_properties)

def geodataframe2rhealpix(gdf, resolution, compact=False, include_properties=True):
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
    return geometry2rhealpix(geometries, resolution, properties_list, compact, include_properties)

def vector2rhealpix(data, resolution, compact=False, output_format='geojson', output_path=None, include_properties=True, **kwargs):
    if not isinstance(resolution, int) or resolution < 0 or resolution > 15:
        raise ValueError(f"Resolution must be in range [0..15], got {resolution}")
    if hasattr(data, 'geometry') and hasattr(data, 'columns'):
        result = geodataframe2rhealpix(data, resolution, compact, include_properties)
    elif isinstance(data, pd.DataFrame):
        result = dataframe2rhealpix(data, resolution, compact, include_properties)
    elif hasattr(data, 'geom_type') or (isinstance(data, list) and len(data) > 0 and hasattr(data[0], 'geom_type')):
        result = geometry2rhealpix(data, resolution, None, compact, include_properties)
    elif isinstance(data, dict) and 'type' in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data['features'])
            result = geodataframe2rhealpix(gdf, resolution, compact, include_properties)
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2rhealpix(gdf, resolution, compact, include_properties)
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
            gdf.to_file("vector2rhealpix.gpkg", driver='GPKG')
            return "vector2rhealpix.gpkg"
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
            gdf.to_file("vector2rhealpix.shp", driver="ESRI Shapefile")
            return "vector2rhealpix.shp"
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile")

def vector2rhealpix_cli():
    parser = argparse.ArgumentParser(description='Convert vector data to rHEALPix grid cells')
    parser.add_argument('-i', '--input', help='Input file path, URL')
    parser.add_argument('-r', '--resolution', type=int, choices=range(0, 16), metavar='[0-15]', help='rHEALPix resolution [0..15] (0=coarsest, 15=finest)')
    parser.add_argument('-c', '--compact', action='store_true', help='Enable rHEALPix compact mode for polygons')
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
        output_path = f"vector2rhealpix{extensions.get(args.format, '')}"
    try:
        result = vector2rhealpix(
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
    vector2rhealpix_cli() 
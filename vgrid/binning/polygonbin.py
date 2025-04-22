import argparse, os, json
from shapely.geometry import shape, Point
from collections import defaultdict, Counter
import statistics
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value

def polygon_bin(polygon_features, point_features, stat, category=None, field_name=None):
    """
    Bins points into arbitrary polygon features and computes statistics.

    :param polygon_features: list of polygon GeoJSON features
    :param point_features: list of point GeoJSON features
    :param stat: statistic type (count, sum, mean, etc.)
    :param category: optional grouping field in properties
    :param field_name: field to use for statistical calculation (except for count)
    :return: list of polygon features with stats in properties
    """
    # Preprocess polygons
    polygons = []
    for feature in polygon_features:
        poly = shape(feature['geometry'])
        if not poly.is_valid:
            continue
        polygons.append((feature, poly))

    # Initialize stat structure per polygon ID
    polygon_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

    for point_feature in tqdm(point_features, desc="Binning points into polygons"):
        point_geom = shape(point_feature['geometry'])
        props = point_feature.get('properties', {})

        if not isinstance(point_geom, Point):
            continue

        for poly_feature, poly_shape in polygons:
            if poly_shape.contains(point_geom):
                poly_id = id(poly_feature)  # Use object ID as unique key
                append_stats_value(polygon_bins, poly_id, props, stat, category, field_name)
                break  # one point belongs to only one polygon

    # Attach stats to polygon properties
    result_features = []
    for poly_feature, poly_shape in polygons:
        poly_id = id(poly_feature)
        categories = polygon_bins.get(poly_id, {})
        out_feature = poly_feature.copy()
        out_feature['properties'] = out_feature.get('properties', {}).copy()

        for cat, values in categories.items():
            key_prefix = '' if category is None else f'{cat}_'

            if stat == 'count':
                out_feature['properties'][f'{key_prefix}count'] = values['count']
            elif stat == 'sum':
                out_feature['properties'][f'{key_prefix}sum'] = sum(values['sum'])
            elif stat == 'min':
                out_feature['properties'][f'{key_prefix}min'] = min(values['min'])
            elif stat == 'max':
                out_feature['properties'][f'{key_prefix}max'] = max(values['max'])
            elif stat == 'mean':
                out_feature['properties'][f'{key_prefix}mean'] = statistics.mean(values['mean'])
            elif stat == 'median':
                out_feature['properties'][f'{key_prefix}median'] = statistics.median(values['median'])
            elif stat == 'std':
                out_feature['properties'][f'{key_prefix}std'] = statistics.stdev(values['std']) if len(values['std']) > 1 else 0
            elif stat == 'var':
                out_feature['properties'][f'{key_prefix}var'] = statistics.variance(values['var']) if len(values['var']) > 1 else 0
            elif stat == 'range':
                out_feature['properties'][f'{key_prefix}range'] = max(values['range']) - min(values['range']) if values['range'] else 0
            elif stat == 'minority':
                freq = Counter(values['values'])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                out_feature['properties'][f'{key_prefix}minority'] = min_item
            elif stat == 'majority':
                freq = Counter(values['values'])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                out_feature['properties'][f'{key_prefix}majority'] = max_item
            elif stat == 'variety':
                out_feature['properties'][f'{key_prefix}variety'] = len(set(values['values']))

        result_features.append(out_feature)

    return result_features

def main(): 
    parser = argparse.ArgumentParser(description="Bin points into polygons and compute statistics")
    parser.add_argument('-point', '--point', type=str, required=True, help="GeoJSON file path (Point or MultiPoint)")
    parser.add_argument('-polygon', '--polygon', type=str, required=True, help="Polygon GeoJSON file path")
    parser.add_argument(
        '-stats', '--statistic', choices=[
            'count', 'min', 'max', 'sum', 'mean', 'median',
            'std', 'var', 'range', 'minority', 'majority', 'variety'
        ], required=True,
        help="Statistic option: choose from count, min, max, sum, mean, median, std, var, range, minority, majority, variety"
    )
    parser.add_argument('-category', '--category', required=False, help="Optional category field for grouping")
    parser.add_argument('-field', '--field', required=False, help="Field name for numeric values")

    args = parser.parse_args()
    
    point_path = args.point
    polygon_path = args.polygon
    stats = args.statistic
    category = args.category
    field_name = args.field

    if not os.path.exists(point_path):
        print(f"Error: The file {point_path} does not exist.")
        return

    if not os.path.exists(polygon_path):
        print(f"Error: The file {polygon_path} does not exist.")
        return

    if stats != 'count' and not field_name:
        print("Error: A field name is required for statistics other than 'count'.")
        return

    with open(point_path, 'r', encoding='utf-8') as f:
        point_data = json.load(f)
    point_features = point_data['features']

    with open(polygon_path, 'r', encoding='utf-8') as f:
        polygon_data = json.load(f)
    polygon_features = polygon_data['features']

    result_features = polygon_bin(polygon_features, point_features, stats, category, field_name)

    out_name = os.path.splitext(os.path.basename(point_path))[0]
    polygon_name = os.path.splitext(os.path.basename(polygon_path))[0]
    out_path = f"{out_name}_bin_{polygon_name}_{stats}.geojson"

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump({"type": "FeatureCollection", "features": result_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")

if __name__ == "__main__":
    main()

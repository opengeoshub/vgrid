import argparse
import os
import json
from collections import defaultdict,Counter
import h3
import statistics
from shapely.geometry import Point, Polygon,shape,mapping
from tqdm import tqdm
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import  geodesic_dggs_to_feature

def _get_default_stat_structure():
    return {
        'count': 0,
        'sum': [],
        'mean': [],
        'min': [],
        'max': [],
        'median': [],
        'std': [],
        'var': [],
        'range': [],
        'values': []  # For variety, minority, majority
    }

def _append_stat_value(h3_bins, h3_id, props, stats, category, field_name=None):
    category_value = props.get(category, "all") if category else "all"
    if h3_id not in h3_bins:
        h3_bins[h3_id] = defaultdict(_get_default_stat_structure)

    if stats == 'count':
        h3_bins[h3_id][category_value]['count'] += 1
    elif stats in ['minority', 'majority', 'variety']:
        value = props.get(field_name or category)
        if value is not None:
            h3_bins[h3_id][category_value]['values'].append(value)
    elif field_name:
        value = props.get(field_name)
        if value is not None:
            try:
                val = float(value)
                h3_bins[h3_id][category_value][stats].append(val)
            except ValueError:
                pass

def point_binning(resolution, point_features, stats, category, field_name, h3_geojson=None):
    h3_bins = defaultdict(lambda: defaultdict(_get_default_stat_structure))
    h3_geometries = {}  # h3_id -> shapely Polygon
    h3_features_map = {}  # h3_id -> original feature (for geometry reuse)

    # 1. Parse h3_geojson features and collect h3_ids and geometries
    if h3_geojson:
        for feature in h3_geojson['features']:
            h3_id = feature['properties'].get('h3')
            if not h3_id:
                continue
            geom = shape(feature['geometry'])
            h3_geometries[h3_id] = geom
            h3_features_map[h3_id] = feature

    # 2. Bin input points to h3_ids
    for feature in tqdm(point_features, desc="Binning points"):
        geom = feature['geometry']
        props = feature.get('properties', {})

        coords_list = []
        if geom['type'] == 'Point':
            coords_list = [geom['coordinates']]
        elif geom['type'] == 'MultiPoint':
            coords_list = geom['coordinates']

        for coords in coords_list:
            point = Point(coords)
            h3_id = h3.latlng_to_cell(point.y, point.x, resolution)

            # Only bin to H3s from the provided GeoJSON
            if h3_geojson and h3_id not in h3_geometries:
                continue

            _append_stat_value(h3_bins, h3_id, props, stats, category, field_name)

    # 3. If no geojson provided, generate geometries for used h3_ids
    if not h3_geojson:
        for h3_id in h3_bins.keys():
            coords = [(lng, lat) for lat, lng in h3.cell_to_boundary(h3_id)]
            geom = Polygon(coords)
            h3_geometries[h3_id] = geom


    # 4. Assemble result features
    result_features = []

    for h3_id, geom in h3_geometries.items():
        props = {}
        categories = h3_bins.get(h3_id, {})

        if not categories:
            # Still include the feature with just the original properties
            feature = {
                "type": "Feature",
                "geometry": mapping(geom),
                "properties": {"h3": h3_id}
            }
            result_features.append(feature)
            continue

        for cat, values in categories.items():
            key_prefix = '' if category is None else f'{cat}_'

            if stats == 'count':
                props[f'{key_prefix}count'] = values['count']
            elif stats == 'sum':
                props[f'{key_prefix}sum'] = sum(values['sum'])
            elif stats == 'mean':
                props[f'{key_prefix}mean'] = statistics.mean(values['mean'])
            elif stats == 'min':
                props[f'{key_prefix}min'] = min(values['min'])
            elif stats == 'max':
                props[f'{key_prefix}max'] = max(values['max'])
            elif stats == 'median':
                props[f'{key_prefix}median'] = statistics.median(values['median'])
            elif stats == 'std':
                props[f'{key_prefix}std'] = statistics.stdev(values['std']) if len(values['std']) > 1 else 0
            elif stats == 'var':
                props[f'{key_prefix}var'] = statistics.variance(values['var']) if len(values['var']) > 1 else 0
            elif stats == 'range':
                props[f'{key_prefix}range'] = max(values['range']) - min(values['range']) if values['range'] else 0
            elif stats == 'minority':
                freq = Counter(values['values'])
                props[f'{key_prefix}minority'] = min(freq.items(), key=lambda x: x[1])[0] if freq else None
            elif stats == 'majority':
                freq = Counter(values['values'])
                props[f'{key_prefix}majority'] = max(freq.items(), key=lambda x: x[1])[0] if freq else None
            elif stats == 'variety':
                props[f'{key_prefix}variety'] = len(set(values['values']))

        feature = {
            "type": "Feature",
            "geometry": mapping(geom),
            "properties": {"h3": h3_id, **props}
        }

        result_features.append(feature)

    return result_features


def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to H3 Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..15]")
    parser.add_argument('-point', '--point', type=str, required=True, help="GeoJSON file path (Point or MultiPoint)")
    parser.add_argument(
        '-stats', '--statistic', choices=[
        'count', 'min', 'max', 'sum', 'mean', 'median',
        'std', 'var', 'range', 'minority', 'majority', 'variety'
        ], required=True,
        help="Statistic option: choose from count, min, max, sum, mean, median,'std', 'var', 'range', 'minority', 'majority', 'variety'"
    )
    parser.add_argument('-category', '--category', required=False, help="Optional category field for grouping")
    parser.add_argument('-field', '--field', required=False, help="Field name for numeric values")
    parser.add_argument('-h3', '--h3', type=str, required=False, help="H3 GeoJSON file with h3 field")

    args = parser.parse_args()

    resolution = args.resolution
    point = args.point
    stats = args.statistic
    category = args.category
    field_name = args.field
    h3_path = args.h3

    if resolution < 0 or resolution > 15:
        print("Error: Please select a resolution in [0..15].")
        return

    if not os.path.exists(point):
        print(f"Error: The file {point} does not exist.")
        return

    if stats != 'count' and not field_name:
        print("Error: A field name is required for statistics other than 'count'.")
        return

    with open(point, 'r', encoding='utf-8') as f:
        point_data = json.load(f)
    point_features = point_data['features']

    h3_geojson = None
    if h3_path:
        if not os.path.exists(h3_path):
            print(f"Error: The H3 file {h3_path} does not exist.")
            return
        with open(h3_path, 'r', encoding='utf-8') as f:
            h3_geojson = json.load(f)

    result_features = point_binning(resolution, point_features, stats, category, field_name, h3_geojson)

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_binning_h3_{resolution}_{stats}.geojson"

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump({"type": "FeatureCollection", "features": result_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")

if __name__ == "__main__":
    main()


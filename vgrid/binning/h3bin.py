import argparse, os, json, statistics
from collections import defaultdict,Counter
import h3
from shapely.geometry import Point, Polygon
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_tructure, append_stats_value
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import  geodesic_dggs_to_feature

def h3_bin(point_features, resolution, stats, category, field_name):
    h3_bins = defaultdict(lambda: defaultdict(get_default_stats_tructure))

    for feature in tqdm(point_features, desc="Binning points"):
        geom = feature['geometry']
        props = feature.get('properties', {})

        if geom['type'] == 'Point':
            point = Point(geom['coordinates'])
            h3_id = h3.latlng_to_cell(point.y, point.x, resolution)
            append_stats_value(h3_bins, h3_id, props, stats, category, field_name)

        elif geom['type'] == 'MultiPoint':
            for coords in geom['coordinates']:
                point = Point(coords)
                h3_id = h3.latlng_to_cell(point.y, point.x, resolution)
                append_stats_value(h3_bins, h3_id, props, stats, category, field_name)

    h3_features = []
    for h3_id, categories in h3_bins.items():
        cell_boundary = h3.cell_to_boundary(h3_id)
        polygon = Polygon([(lon, lat) for lat, lon in fix_h3_antimeridian_cells(cell_boundary)])

        if not polygon.is_valid:
            continue

        h3_feature = geodesic_dggs_to_feature("h3", h3_id, resolution, polygon, 5 if h3.is_pentagon(h3_id) else 6)

        for cat, values in categories.items():
            key_prefix = '' if category is None else f'{cat}_'

            if stats == 'count':
                h3_feature['properties'][f'{key_prefix}count'] = values['count']
            elif stats == 'sum':
                h3_feature['properties'][f'{key_prefix}sum'] = sum(values['sum'])
            elif stats == 'min':
                h3_feature['properties'][f'{key_prefix}min'] = min(values['min'])
            elif stats == 'max':
                h3_feature['properties'][f'{key_prefix}max'] = max(values['max'])
            elif stats == 'mean':
                h3_feature['properties'][f'{key_prefix}mean'] = statistics.mean(values['mean'])
            elif stats == 'median':
                h3_feature['properties'][f'{key_prefix}median'] = statistics.median(values['median'])
            elif stats == 'std':
                h3_feature['properties'][f'{key_prefix}std'] = statistics.stdev(values['std']) if len(values['std']) > 1 else 0
            elif stats == 'var':
                h3_feature['properties'][f'{key_prefix}var'] = statistics.variance(values['var']) if len(values['var']) > 1 else 0
            elif stats == 'range':
                h3_feature['properties'][f'{key_prefix}range'] = max(values['range']) - min(values['range']) if values['range'] else 0
            
            elif stats == 'minority':
                freq = Counter(values['values'])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                h3_feature['properties'][f'{key_prefix}minority'] = min_item

            elif stats == 'majority':
                freq = Counter(values['values'])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                h3_feature['properties'][f'{key_prefix}majority'] = max_item

            elif stats == 'variety':
                h3_feature['properties'][f'{key_prefix}variety'] = len(set(values['values']))


        h3_features.append(h3_feature)

    return h3_features


def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to H3 Grid")
    parser.add_argument('-point', '--point', type=str, required=True, help="GeoJSON file path (Point or MultiPoint)")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..15]")
    parser.add_argument(
            '-stats', '--statistics',
            choices=[
                'count', 'min', 'max', 'sum', 'mean', 'median',
                'std', 'var', 'range', 'minority', 'majority', 'variety'
            ],
            required=True,
            help="Statistic option"
        )

    parser.add_argument('-category', '--category', required=False, help="Optional category field for grouping")
    parser.add_argument('-field', '--field', required=False, help="Numeric field to compute statistics")

    args = parser.parse_args()

    resolution = args.resolution
    point = args.point
    stats = args.statistics
    category = args.category
    field_name = args.field

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
    h3_features = h3_bin(point_features,resolution, stats, category, field_name)

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_bin_h3_{resolution}_{stats}.geojson"

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump({"type": "FeatureCollection", "features": h3_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()
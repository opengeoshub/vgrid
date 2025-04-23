import argparse, os, json, statistics
from collections import defaultdict,Counter
from shapely.geometry import Point
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value
from vgrid.generator.settings import  geodesic_dggs_to_feature
from vgrid.utils import qtm

def qtm_bin(point_features, resolution, stats, category, field_name):
    qtm_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

    for feature in tqdm(point_features, desc="Binning points"):
        geom = feature['geometry']
        props = feature.get('properties', {})

        if geom['type'] == 'Point':
            point = Point(geom['coordinates'])
            qtm_id = qtm.latlon_to_qtm_id(point.y, point.x, resolution)
            append_stats_value(qtm_bins, qtm_id, props, stats, category, field_name)

        elif geom['type'] == 'MultiPoint':
            for coords in geom['coordinates']:
                point = Point(coords)
                qtm_id = qtm.latlon_to_qtm_id(point.y, point.x, resolution)
                append_stats_value(qtm_bins, qtm_id, props, stats, category, field_name)

    qtm_features = []
    for qtm_id, categories in qtm_bins.items():
        facet = qtm.qtm_id_to_facet(qtm_id)
        cell_polygon = qtm.constructGeometry(facet)    

        if not cell_polygon.is_valid:
            continue
        
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature("qtm", qtm_id, resolution, cell_polygon,num_edges)

        for cat, values in categories.items():
            key_prefix = '' if category is None else f'{cat}_'

            if stats == 'count':
                qtm_feature['properties'][f'{key_prefix}count'] = values['count']
            elif stats == 'sum':
                qtm_feature['properties'][f'{key_prefix}sum'] = sum(values['sum'])
            elif stats == 'min':
                qtm_feature['properties'][f'{key_prefix}min'] = min(values['min'])
            elif stats == 'max':
                qtm_feature['properties'][f'{key_prefix}max'] = max(values['max'])
            elif stats == 'mean':
                qtm_feature['properties'][f'{key_prefix}mean'] = statistics.mean(values['mean'])
            elif stats == 'median':
                qtm_feature['properties'][f'{key_prefix}median'] = statistics.median(values['median'])
            elif stats == 'std':
                qtm_feature['properties'][f'{key_prefix}std'] = statistics.stdev(values['std']) if len(values['std']) > 1 else 0
            elif stats == 'var':
                qtm_feature['properties'][f'{key_prefix}var'] = statistics.variance(values['var']) if len(values['var']) > 1 else 0
            elif stats == 'range':
                qtm_feature['properties'][f'{key_prefix}range'] = max(values['range']) - min(values['range']) if values['range'] else 0
            
            elif stats == 'minority':
                freq = Counter(values['values'])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                qtm_feature['properties'][f'{key_prefix}minority'] = min_item

            elif stats == 'majority':
                freq = Counter(values['values'])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                qtm_feature['properties'][f'{key_prefix}majority'] = max_item

            elif stats == 'variety':
                qtm_feature['properties'][f'{key_prefix}variety'] = len(set(values['values']))


        qtm_features.append(qtm_feature)

    return qtm_features


def main():
    parser = argparse.ArgumentParser(description="Binning point to QTM DGGS")
    parser.add_argument('-point', '--point', type=str, required=True, help="GeoJSON file path (Point or MultiPoint)")
    parser.add_argument('-r', '--resolution', type=int, default=14, help="Resolution of the grid [1..24]")
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

    if resolution < 1 or resolution > 24:
        print("Error: Please select a resolution in [1..24].")
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
    qtm_features = qtm_bin(point_features,resolution, stats, category, field_name)

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_bin_qtm_{resolution}_{stats}.geojson"

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump({"type": "FeatureCollection", "features": qtm_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()
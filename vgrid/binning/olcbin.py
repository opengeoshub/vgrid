import argparse, os, json, statistics
from collections import defaultdict,Counter
from shapely.geometry import Point,Polygon
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value
from vgrid.generator.settings import  graticule_dggs_to_feature
from vgrid.conversion.latlon2dggs import latlon2olc
from vgrid.utils import olc

def olc_bin(point_features, resolution, stats, category, field_name):
    olc_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

    for feature in tqdm(point_features, desc="Binning points"):
        geom = feature['geometry']
        props = feature.get('properties', {})

        if geom['type'] == 'Point':
            point = Point(geom['coordinates'])
            olc_id = latlon2olc(point.y, point.x, resolution)
            append_stats_value(olc_bins, olc_id, props, stats, category, field_name)

        elif geom['type'] == 'MultiPoint':
            for coords in geom['coordinates']:
                point = Point(coords)
                olc_id = latlon2olc(point.y, point.x, resolution)
                append_stats_value(olc_bins, olc_id, props, stats, category, field_name)

    olc_features = []
    for olc_id, categories in olc_bins.items():
        coord = olc.decode(olc_id)    
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi        
        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])  

        if not cell_polygon.is_valid:
            continue
        
        olc_feature = graticule_dggs_to_feature("olc", olc_id, resolution, cell_polygon)

        for cat, values in categories.items():
            key_prefix = '' if category is None else f'{cat}_'

            if stats == 'count':
                olc_feature['properties'][f'{key_prefix}count'] = values['count']
            elif stats == 'sum':
                olc_feature['properties'][f'{key_prefix}sum'] = sum(values['sum'])
            elif stats == 'min':
                olc_feature['properties'][f'{key_prefix}min'] = min(values['min'])
            elif stats == 'max':
                olc_feature['properties'][f'{key_prefix}max'] = max(values['max'])
            elif stats == 'mean':
                olc_feature['properties'][f'{key_prefix}mean'] = statistics.mean(values['mean'])
            elif stats == 'median':
                olc_feature['properties'][f'{key_prefix}median'] = statistics.median(values['median'])
            elif stats == 'std':
                olc_feature['properties'][f'{key_prefix}std'] = statistics.stdev(values['std']) if len(values['std']) > 1 else 0
            elif stats == 'var':
                olc_feature['properties'][f'{key_prefix}var'] = statistics.variance(values['var']) if len(values['var']) > 1 else 0
            elif stats == 'range':
                olc_feature['properties'][f'{key_prefix}range'] = max(values['range']) - min(values['range']) if values['range'] else 0
            
            elif stats == 'minority':
                freq = Counter(values['values'])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                olc_feature['properties'][f'{key_prefix}minority'] = min_item

            elif stats == 'majority':
                freq = Counter(values['values'])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                olc_feature['properties'][f'{key_prefix}majority'] = max_item

            elif stats == 'variety':
                olc_feature['properties'][f'{key_prefix}variety'] = len(set(values['values']))


        olc_features.append(olc_feature)

    return olc_features


def main():
    parser = argparse.ArgumentParser(description="Binning point to OLC DGGS")
    parser.add_argument('-point', '--point', type=str, required=True, help="GeoJSON file path (Point or MultiPoint)")
    parser.add_argument(
            '-r', '--resolution',
            type=int,
            choices=[2, 4, 6, 8, 10, 11, 12, 13, 14, 15],
            default=8,
            help="Resolution of the OLC DGGS (choose from 2, 4, 6, 8, 10, 11, 12, 13, 14, 15)"
        )
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

    if not os.path.exists(point):
        print(f"Error: The file {point} does not exist.")
        return

    if stats != 'count' and not field_name:
        print("Error: A field name is required for statistics other than 'count'.")
        return

    with open(point, 'r', encoding='utf-8') as f:
        point_data = json.load(f)

    point_features = point_data['features']
    olc_features = olc_bin(point_features,resolution, stats, category, field_name)

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_bin_olc_{resolution}_{stats}.geojson"

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump({"type": "FeatureCollection", "features": olc_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()
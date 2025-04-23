import argparse, os, re, json, statistics
from collections import defaultdict,Counter
from shapely.geometry import Point,Polygon
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value
from vgrid.generator.settings import  graticule_dggs_to_feature
from vgrid.conversion.latlon2dggs import latlon2tilecode
from vgrid.utils import mercantile

def tilecode_bin(point_features, resolution, stats, category, field_name):
    tilecode_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

    for feature in tqdm(point_features, desc="Binning points"):
        geom = feature['geometry']
        props = feature.get('properties', {})

        if geom['type'] == 'Point':
            point = Point(geom['coordinates'])
            tilecode_id = latlon2tilecode(point.y, point.x, resolution)
            append_stats_value(tilecode_bins, tilecode_id, props, stats, category, field_name)

        elif geom['type'] == 'MultiPoint':
            for coords in geom['coordinates']:
                point = Point(coords)
                tilecode_id = latlon2tilecode(point.y, point.x, resolution)
                append_stats_value(tilecode_bins, tilecode_id, props, stats, category, field_name)

    tilecode_features = []
    for tilecode_id, categories in tilecode_bins.items():
        match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)    
        # Convert matched groups to integers
        z = int(match.group(1))
        x = int(match.group(2))
        y = int(match.group(3))
        bounds = mercantile.bounds(x, y, z)    
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])

        if not cell_polygon.is_valid:
            continue
        
        tilecode_feature = graticule_dggs_to_feature("tilecode", tilecode_id, resolution, cell_polygon)

        for cat, values in categories.items():
            key_prefix = '' if category is None else f'{cat}_'

            if stats == 'count':
                tilecode_feature['properties'][f'{key_prefix}count'] = values['count']
            elif stats == 'sum':
                tilecode_feature['properties'][f'{key_prefix}sum'] = sum(values['sum'])
            elif stats == 'min':
                tilecode_feature['properties'][f'{key_prefix}min'] = min(values['min'])
            elif stats == 'max':
                tilecode_feature['properties'][f'{key_prefix}max'] = max(values['max'])
            elif stats == 'mean':
                tilecode_feature['properties'][f'{key_prefix}mean'] = statistics.mean(values['mean'])
            elif stats == 'median':
                tilecode_feature['properties'][f'{key_prefix}median'] = statistics.median(values['median'])
            elif stats == 'std':
                tilecode_feature['properties'][f'{key_prefix}std'] = statistics.stdev(values['std']) if len(values['std']) > 1 else 0
            elif stats == 'var':
                tilecode_feature['properties'][f'{key_prefix}var'] = statistics.variance(values['var']) if len(values['var']) > 1 else 0
            elif stats == 'range':
                tilecode_feature['properties'][f'{key_prefix}range'] = max(values['range']) - min(values['range']) if values['range'] else 0
            
            elif stats == 'minority':
                freq = Counter(values['values'])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                tilecode_feature['properties'][f'{key_prefix}minority'] = min_item

            elif stats == 'majority':
                freq = Counter(values['values'])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                tilecode_feature['properties'][f'{key_prefix}majority'] = max_item

            elif stats == 'variety':
                tilecode_feature['properties'][f'{key_prefix}variety'] = len(set(values['values']))


        tilecode_features.append(tilecode_feature)

    return tilecode_features


def main():
    parser = argparse.ArgumentParser(description="Binning point to tilecode DGGS")
    parser.add_argument('-point', '--point', type=str, required=True, help="GeoJSON file path (Point or MultiPoint)")
    parser.add_argument('-r', '--resolution', type=int, default=15, help="Resolution of the grid [0..29]")
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

    if resolution < 0 or resolution > 29:
        print("Error: Please select a resolution in [0..29].")
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
    tilecode_features = tilecode_bin(point_features,resolution, stats, category, field_name)

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_bin_tilecode_{resolution}_{stats}.geojson"

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump({"type": "FeatureCollection", "features": tilecode_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()
#Reference: https://github.com/google/open-location-code/tree/main/tile_server/gridserver
# https://github.com/google/open-location-code

import json
import argparse
import tiles_util.utils.geocode.olc as olc
from tqdm import tqdm
from shapely.geometry import Polygon

def generate_all_olcs(length):
    olc_chars = '23456789CFGHJMPQRVWX'
    if length < 2:
        raise ValueError("OLC length should be at least 2.")

    def olc_generator(prefix, depth):
        if depth == length:
            yield prefix
        else:
            for char in olc_chars:
                yield from olc_generator(prefix + char, depth + 1)

    return olc_generator("", 0)

def create_geojson_for_olc(olc_code, length):
    decoded = olc.decode(olc_code)
    coordinates = [
        [decoded.longitudeLo, decoded.latitudeLo],
        [decoded.longitudeLo, decoded.latitudeHi],
        [decoded.longitudeHi, decoded.latitudeHi],
        [decoded.longitudeHi, decoded.latitudeLo],
        [decoded.longitudeLo, decoded.latitudeLo]
    ]
    polygon = Polygon(coordinates)
    polygon_feature = {
        "type": "Feature",
        "geometry": {
            "type": "Polygon",
            "coordinates": [coordinates]
        },
        "properties": {
            "name": f"polygon_{length}",
            "pluscode": olc_code
        }
    }  
    return polygon_feature

def is_within_bounding_box(decoded, bbox):
    return (decoded.longitudeLo < bbox[2] and decoded.longitudeHi > bbox[0] and
            decoded.latitudeLo < bbox[3] and decoded.latitudeHi > bbox[1])

def generate_geojson_for_olc_length(length, bbox):
    features = []
    total_codes = 20 ** length  # Total number of possible codes of the given length
    for olc_code in tqdm(generate_all_olcs(length), total=total_codes, desc="Generating GeoJSON"):
        decoded = olc.decode(olc_code)
        if is_within_bounding_box(decoded, bbox):
            polygon_feature = create_geojson_for_olc(olc_code, length)
            features.append(polygon_feature)   
    geojson = {
        "type": "FeatureCollection",
        "features": features
    }
    
    return geojson

def main():
    parser = argparse.ArgumentParser(description="Generate GeoJSON with OLC codes and centroids.")
    parser.add_argument('-l', '--length', type=int, help="Length of the plus code [2,4,8]") 
    # length = 2: zoomlevel 0 - 5
    # length = 4: zoomlevel 6 - 10
    # length = 8: zoomlevel 11 - >11
    args = parser.parse_args()
    
    length = args.length
    bbox = [-180, -85.051129, 180, 85.051129]
    
    geojson = generate_geojson_for_olc_length(length, bbox)
    
    filename = f'pluscode_{length}.geojson'
    with open(filename, 'w') as f:
        json.dump(geojson, f, indent=2)
    print(f"GeoJSON file saved as {filename}")

if __name__ == "__main__":
    main()
import argparse
import re
import json
from shapely.geometry import mapping,Polygon
from tqdm import tqdm
from vgrid.utils import mercantile
from pyproj import Geod
geod = Geod(ellps="WGS84")
import locale
current_locale = locale.getlocale()  # Get the current locale setting
locale.setlocale(locale.LC_ALL,current_locale)  # Use the system's default locale

max_cells = 1_000_000

def generate_grid(resolution,bbox=None):
    features = []
    min_lon, min_lat, max_lon, max_lat = bbox or [-180.0, -85.05112878,180.0,85.05112878]  
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    for tile in tqdm(tiles, desc=f"Processing tiles at zoom level {resolution}:", unit=" cells"):
        z, x, y = tile.z, tile.x, tile.y
        tilecode = f"z{tile.z}x{tile.x}y{tile.y}"
        bounds = mercantile.bounds(x, y, z)
        if bounds:
            # Create the bounding box coordinates for the polygon
            min_lat, min_lon = bounds.south, bounds.west
            max_lat, max_lon = bounds.north, bounds.east
            
            quadkey = mercantile.quadkey(tile)

            center_lat = round((min_lat + max_lat) / 2,7)
            center_lon = round((min_lon + max_lon) / 2,7)
            
            cell_polygon = Polygon([
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat]   # Closing the polygon (same as the first point)
            ])
            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)  # Area in square meters     
            # Calculate width (longitude difference at a constant latitude)
            cell_width = round(geod.line_length([min_lon, max_lon], [min_lat, min_lat]),3)
            # Calculate height (latitude difference at a constant longitude)
            cell_height = round(geod.line_length([min_lon, min_lon], [min_lat, max_lat]),3)

            feature = {
                "type": "Feature",
                "geometry": mapping(cell_polygon),          
                "properties": {
                    "tilecode": tilecode,  # Include the OLC as a property
                    "quadkey": quadkey,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_area": cell_area,
                    "cell_width": cell_width,
                    "cell_height": cell_height,
                    "resolution": z  # Using the code length as precision
                }
            }            

        features.append(feature)

    geojson_features = {
        'type': 'FeatureCollection',
        'features': features
    }

    return geojson_features
        
def main():
    parser = argparse.ArgumentParser(description='Create a debug grid representing the XYZ vector tile scheme as a GeoJSON file.')
    parser.add_argument('-r', '--resolution', type=int, required=True, help='zoom level/ resolution= [0..26]')
    parser.add_argument('-b', '--bbox', type=float, nargs=4,  help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)") 

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180.0, -85.05112878,  180.0, 85.05112878]
    
    if resolution < 0 or resolution > 26:
        print(f"Please select a resolution in [0..26] range and try again ")
        return
    
    if bbox == [-180.0, -85.05112878,  180.0, 85.05112878]:  
        num_cells =  4**resolution
        if num_cells > max_cells:
            print(
                f"The selected resolution will generate "
                f"{locale.format_string('%d', num_cells, grouping=True)} cells, "
                f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}."
            )
            print("Please select a smaller resolution and try again.")
            return    
    
    geojson_features = generate_grid(resolution, bbox)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"tile_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")

if __name__ == '__main__':
    main()

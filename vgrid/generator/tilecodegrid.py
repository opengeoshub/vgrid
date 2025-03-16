import argparse
import re
import json
from shapely.geometry import mapping,Polygon
from tqdm import tqdm
from vgrid.utils import mercantile
from vgrid.generator.settings import max_cells, graticule_dggs_to_feature

def generate_grid(resolution,bbox):
    tilecode_features = []
    min_lon, min_lat, max_lon, max_lat = bbox # or [-180.0, -85.05112878,180.0,85.05112878]  
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    
    for tile in tqdm(tiles, desc=f"Processing tiles at zoom level {resolution}", unit=" cells"):
        z, x, y = tile.z, tile.x, tile.y
        tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
        bounds = mercantile.bounds(x, y, z)
        if bounds:
            # Create the bounding box coordinates for the polygon
            min_lat, min_lon = bounds.south, bounds.west
            max_lat, max_lon = bounds.north, bounds.east          
            
            cell_polygon = Polygon([
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat]   # Closing the polygon (same as the first point)
            ])
            tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id,resolution,cell_polygon)   
            tilecode_features.append(tilecode_feature)

    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }

        
def main():
    parser = argparse.ArgumentParser(description='Generate Tilecode grid.')
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
                f"The selected resolution will generate {num_cells} cells "
                f"which exceeds the limit of {max_cells}."
            )
            print("Please select a smaller resolution and try again.")
            return    
    
    geojson_features = generate_grid(resolution, bbox)
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"tilecode_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")

if __name__ == '__main__':
    main()

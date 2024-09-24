import argparse
import os
from vgrid.utils.mapbox_vector_tile import encode
import vgrid.utils.mercantile as mercantile
from shapely.geometry import box, mapping
import sqlite3
import gzip
from tqdm import tqdm

def create_tile(z, x, y):
    try:
        flip_y = (2 ** z - 1) - y
        tile_geometry = box(0, 0, 4096, 4096)
        quadkey = mercantile.quadkey(x, flip_y, z)
        properties = {
            'vcode': f'z{z}x{x}y{flip_y}',
            'tmscode': f'z{z}x{x}y{y}',
            'vname': f'{z}/{x}/{flip_y}',
            'tms_name': f'{z}/{x}/{y}',
            'quadkey': quadkey
        }

        feature = {
            'geometry': mapping(tile_geometry),
            'properties': properties
        }

        tile_data = {
            'name': 'vectortile_grid',
            'features': [feature]
        }

        tile_data_encoded = encode(tile_data)
        tile_data_encoded_gzipped = gzip.compress(tile_data_encoded)
        return tile_data_encoded_gzipped
    except Exception as e:
        print(f"Error creating tile: {e}")
        raise

def create_mbtiles(tiles, output_mbtiles, min_zoom, max_zoom):
    try:
        if os.path.exists(output_mbtiles):
            os.remove(output_mbtiles)

        conn = sqlite3.connect(output_mbtiles)
        cursor = conn.cursor()

        cursor.execute('''CREATE TABLE metadata (name TEXT, value TEXT);''')
        cursor.execute('''CREATE UNIQUE INDEX name ON metadata (name);''')
        cursor.execute('''CREATE TABLE tiles (zoom_level INTEGER, tile_column INTEGER, tile_row INTEGER, tile_data BLOB);''')
        cursor.execute('''CREATE UNIQUE INDEX tile_index ON tiles(zoom_level, tile_column, tile_row);''')

        cursor.execute("INSERT INTO metadata (name, value) VALUES ('name', 'Vectortile Grid');")
        cursor.execute("INSERT INTO metadata (name, value) VALUES ('description', 'A vector tile grid created by tiles_util.vectortilegrid');")
        cursor.execute("INSERT INTO metadata (name, value) VALUES ('type', 'overlay');")
        cursor.execute("INSERT INTO metadata (name, value) VALUES ('version', '1');")
        cursor.execute("INSERT INTO metadata (name, value) VALUES ('format', 'pbf');")
        cursor.execute(f"INSERT INTO metadata (name, value) VALUES ('minzoom', '{min_zoom}');")
        cursor.execute(f"INSERT INTO metadata (name, value) VALUES ('maxzoom', '{max_zoom}');")
        cursor.execute("INSERT INTO metadata (name, value) VALUES ('bounds', '-180.000000,-85.051129,180.000000,85.051129');")

        for tile in tqdm(tiles, desc="Creating tiles"):
            z, x, y = tile.z, tile.x, tile.y
            flip_y = (2 ** z - 1) - y
            tile_data = create_tile(z, x, y)
            cursor.execute('INSERT INTO tiles (zoom_level, tile_column, tile_row, tile_data) VALUES (?, ?, ?, ?);', 
                           (z, x, y, tile_data))

        conn.commit()
    except Exception as e:
        print(f"Error creating MBTiles: {e}")
    finally:
        conn.close()

def main():
    parser = argparse.ArgumentParser(description='Create a debug grid representing the XYZ vector tile scheme as an MBTiles file.')
    parser.add_argument('-o', '--output', required=True, help='Output MBTiles file')
    parser.add_argument('-minzoom','--minzoom', type=int, required=True, help='Minimum zoom level')
    parser.add_argument('-maxzoom', '--maxzoom',type=int, required=True, help='Maximum zoom level')
    args = parser.parse_args()

    if args.minzoom < 0 or args.maxzoom < args.minzoom:
        raise ValueError("minzoom must be non-negative and maxzoom must be greater than or equal to minzoom")

    min_latitude = -85.05112878
    max_latitude = 85.05112878
    min_longitude = -180.0
    max_longitude = 180.0

    all_tiles = []

    for zoom_level in range(args.minzoom, args.maxzoom + 1):
        list_tiles = list(mercantile.tiles(min_longitude, min_latitude, max_longitude, max_latitude, zoom_level))
        all_tiles.extend(list_tiles)

    create_mbtiles(all_tiles, args.output, args.minzoom, args.maxzoom)

if __name__ == '__main__':
    main()

from . import mercantile
import math
from shapely.geometry import Polygon
from shapely.ops import transform
import pyproj
from shapely.geometry import box
import argparse
import string
import re

def haversine(lat1, lon1, lat2, lon2):
    # Radius of the Earth in meters
    R = 6371000  

    # Convert latitude and longitude from degrees to radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(delta_phi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return R * c  # Distance in meters


# Define the character set excluding 'z', 'x', and 'y'
characters = string.digits + string.ascii_uppercase + string.ascii_lowercase.replace('z', '').replace('x', '').replace('y', '')
base = len(characters)

def tile_encode(num):
    if num == 0:
        return characters[0]
    
    encoded = []
    while num > 0:
        num, remainder = divmod(num, base)
        encoded.append(characters[remainder])
    
    return ''.join(reversed(encoded))

def tile_decode(encoded):
    num = 0
    for char in encoded:
        num = num * base + characters.index(char)
    return num

def tile_encode_cli():
    parser = argparse.ArgumentParser(description='Encode a number using a custom base encoding (excluding z, x, y).')
    parser.add_argument('number', type=int, help='The number to encode (0-9999999).')
    args = parser.parse_args()
    
    if args.number < 0 or args.number > 9999999:
        print("Error: The number must be between 0 and 9999999.")
        return
    
    encoded_value = tile_encode(args.number)
    print(f"Encoded: {encoded_value}")

def tile_decode_cli():
    parser = argparse.ArgumentParser(description='Decode a custom base encoded string (excluding z, x, y).')
    parser.add_argument('encoded', type=str, help='The encoded string to decode.')
    args = parser.parse_args()
    
    try:
        decoded_value = tile_decode(args.encoded)
        print(f"Decoded: {decoded_value}")
    except ValueError:
        print("Error: The provided encoded string is invalid.")


def tilecode2bbox(tilecode):
    """
    Converts a tilecode (e.g., 'z8x11y14') to a Polygon geometry
    representing the tile's bounds and includes the original tilecode as a property.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        dict: A polygon geometry and tilecode as a property.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile in (west, south, east, north)
    bounds = mercantile.bounds(x, y, z)

    # Create the coordinates of the polygon using the bounds
    polygon_coords = [
        [bounds.west, bounds.south],  # Bottom-left
        [bounds.east, bounds.south],  # Bottom-right
        [bounds.east, bounds.north],  # Top-right
        [bounds.west, bounds.north],  # Top-left
        [bounds.west, bounds.south]   # Closing the polygon
    ]

    return polygon_coords


def zxy2tilecode(z, x, y):
    """
    Converts z, x, and y values to a string formatted as 'zXxYyZ'.

    Args:
        z (int): The zoom level.
        x (int): The x coordinate.
        y (int): The y coordinate.

    Returns:
        str: A string formatted as 'zXxYyZ'.
    """
    return f'z{z}x{x}y{y}'

def tilecode2zxy(tilecode):
    """
    Parses a string formatted as 'zXxYyZ' to extract z, x, and y values.

    Args:
        tilecode (str): A string formatted like 'z8x11y14'.

    Returns:
        tuple: A tuple containing (z, x, y) as integers.
    """
    # Regular expression to capture numbers after z, x, and y
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    
    if match:
        # Extract and convert matched groups to integers
        z = int(match.group(1))
        x = int(match.group(2))
        y = int(match.group(3))
        return z, x, y
    else:
        # Raise an error if the format does not match
        raise ValueError("Invalid format. Expected format: 'zXxYyZ'")

def latlon2tilecode(lat, lon, zoom):
    """
    Converts latitude, longitude, and zoom level to a tile code ('tilecode') of the format 'zXxYyZ'.

    Args:
        lat (float): Latitude of the point.
        lon (float): Longitude of the point.
        zoom (int): Zoom level.

    Returns:
        str: A string representing the tile code in the format 'zXxYyZ'.
    """
    # Get the tile coordinates (x, y) for the given lat, lon, and zoom level
    tile = mercantile.tile(lon, lat, zoom)
    
    # Format the tile coordinates into the tilecode string
    tilecode = f"z{tile.z}x{tile.x}y{tile.y}"
    
    return tilecode

def tilecode2latlon(tilecode):
    """
    Calculates the center latitude and longitude of a tile given its tilecode.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        tuple: A tuple containing the latitude and longitude of the tile's center.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile
    bounds = mercantile.bounds(x, y, z)

    # Calculate the center of the tile
    center_longitude = (bounds.west + bounds.east) / 2
    center_latitude = (bounds.south + bounds.north) / 2

    return [center_latitude, center_longitude]

def tilecode2quadkey(tilecode):
    """
    Converts a tilecode (e.g., 'z23x6668288y3948543') to a quadkey using mercantile.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        str: Quadkey corresponding to the tilecode.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Use mercantile to get the quadkey
    tile = mercantile.Tile(x, y, z)
    quadkey = mercantile.quadkey(tile)

    return quadkey

def quadkey2tilecode(quadkey):
    """
    Converts a quadkey to a tilecode (e.g., 'z23x6668288y3948543') using mercantile.

    Args:
        quadkey (str): The quadkey string.

    Returns:
        str: tilecode in the format 'zXxYyZ'.
    """
    # Decode the quadkey to get the tile coordinates and zoom level
    tile = mercantile.quadkey_to_tile(quadkey)
    
    # Format as tilecode
    tilecode = f"z{tile.z}x{tile.x}y{tile.y}"

    return tilecode

def tilecode_cell_area(tilecode):
    """
    Calculates the area in square meters of a tile given its tilecode.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        float: The area of the tile in square meters.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile
    bounds = mercantile.bounds(x, y, z)

    # Define the polygon from the bounds
    polygon_coords = [
        [bounds.west, bounds.south],  # Bottom-left
        [bounds.east, bounds.south],  # Bottom-right
        [bounds.east, bounds.north],  # Top-right
        [bounds.west, bounds.north],  # Top-left
        [bounds.west, bounds.south]   # Closing the polygon
    ]
    polygon = Polygon(polygon_coords)

    # Project the polygon to a metric CRS (e.g., EPSG:3857) to calculate area in square meters
    project = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True).transform
    metric_polygon = transform(project, polygon)

    # Calculate the area in square meters
    area = metric_polygon.area

    return area

def tilecode_cell_length(tilecode):
    """
    Calculates the length of the edge of a square tile given its tilecode.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        float: The length of the edge of the tile in meters.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile
    bounds = mercantile.bounds(x, y, z)

    # Define the coordinates of the polygon
    polygon_coords = [
        [bounds.west, bounds.south],  # Bottom-left
        [bounds.east, bounds.south],  # Bottom-right
        [bounds.east, bounds.north],  # Top-right
        [bounds.west, bounds.north],  # Top-left
        [bounds.west, bounds.south]   # Closing the polygon
    ]
    polygon = Polygon(polygon_coords)

    # Project the polygon to a metric CRS (e.g., EPSG:3857) to calculate edge length in meters
    project = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True).transform
    metric_polygon = transform(project, polygon)

    # Calculate the length of the edge of the square
    edge_length = metric_polygon.exterior.length / 4  # Divide by 4 for the length of one edge

    return edge_length

def tilecode2tilebound(tilecode):
    """
    Converts a tilecode (e.g., 'z23x6668288y3948543') to its bounding box using mercantile.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        dict: Bounding box with 'west', 'south', 'east', 'north' coordinates.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Use mercantile to get the bounds
    tile = mercantile.Tile(x, y, z)
    bounds = mercantile.bounds(tile)

    # Convert bounds to a dictionary
    bounds_dict = {
        'west': bounds[0],
        'south': bounds[1],
        'east': bounds[2],
        'north': bounds[3]
    }

    return bounds_dict

def tilecode2bound(tilecode):
    """
    Converts a tilecode (e.g., 'z23x6668288y3948543') to its bounding box using mercantile.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        list: Bounding box in the format [left, bottom, right, top].
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Convert tile coordinates to Mercator bounds
    bounds = mercantile.bounds(mercantile.Tile(x, y, z))

    # Return bounds as a list in [left, bottom, right, top] format
    return [bounds[0], bounds[1], bounds[2], bounds[3]]

def tilecode2wktbound(tilecode):
    """
    Converts a tilecode (e.g., 'z23x6668288y3948543') to its bounding box in OGC WKT format using mercantile.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        str: Bounding box in OGC WKT format.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Use mercantile to get the bounds
    tile = mercantile.Tile(x, y, z)
    bounds = mercantile.bounds(tile)

    # Convert bounds to WKT POLYGON format
    wkt = f"POLYGON(({bounds[0]} {bounds[1]}, {bounds[0]} {bounds[3]}, {bounds[2]} {bounds[3]}, {bounds[2]} {bounds[1]}, {bounds[0]} {bounds[1]}))"

    return wkt

def tilecode_list(zoom):
    """
    Lists all tilecodes at a specific zoom level using mercantile.

    Args:
        zoom (int): The zoom level.

    Returns:
        list: A list of tilecodes for the specified zoom level.
    """
    # Get the maximum number of tiles at the given zoom level
    num_tiles = 2 ** zoom

    tilecodes = []
    for x in range(num_tiles):
        for y in range(num_tiles):
            # Create a tile object
            tile = mercantile.Tile(x, y, zoom)
            # Convert tile to tilecode
            tilecode = f"z{tile.z}x{tile.x}y{tile.y}"
            tilecodes.append(tilecode)
    
    return tilecodes


def tilecode_children(tilecode):
    """
    Lists all child tiles of a given tilecode at the next zoom level.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        list: A list of tilecodes representing the four child tiles.
    """
    # Extract z, x, y from the tilecode
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Calculate the next zoom level
    z_next = z + 1

    # Calculate the coordinates of the four child tiles
    children = [
        f"z{z_next}x{2*x}y{2*y}",       # Top-left child
        f"z{z_next}x{2*x+1}y{2*y}",     # Top-right child
        f"z{z_next}x{2*x}y{2*y+1}",     # Bottom-left child
        f"z{z_next}x{2*x+1}y{2*y+1}"    # Bottom-right child
    ]

    return children

def tilecode_parent(tilecode):
    """
    Finds the parent tile of a given tilecode at the current zoom level.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ', where X, Y, and Z are integers.

    Returns:
        str: The tilecode of the parent tile.
    """
    # Extract z, x, y from the tilecode
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Calculate the parent zoom level
    if z == 0:
        raise ValueError("No parent exists for zoom level 0.")

    z_parent = z - 1

    # Calculate the coordinates of the parent tile
    x_parent = x // 2
    y_parent = y // 2

    # Format the parent tile's tilecode
    parent_tilecode = f"z{z_parent}x{x_parent}y{y_parent}"

    return parent_tilecode

def tilecode_siblings(tilecode):
    """
    Lists all sibling tiles of a given tilecode at the same zoom level.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        list: A list of tilecodes representing the sibling tiles, excluding the input tilecode.
    """
    # Extract z, x, y from the tilecode
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Calculate the parent tile's coordinates
    if z == 0:
        # The root tile has no siblings
        return []

    z_parent = z - 1
    x_parent = x // 2
    y_parent = y // 2

    # Get all children of the parent tile
    parent_tilecode = f"z{z_parent}x{x_parent}y{y_parent}"
    children = tilecode_children(parent_tilecode)

    # Exclude the input tilecode from the list of siblings
    siblings = [child for child in children if child != tilecode]

    return siblings

def tilecode_neighbors(tilecode):
    """
    Finds the neighboring tilecodes of a given tilecode.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        list: A list of neighboring tilecodes.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Calculate the neighboring tiles (including the tile itself)
    neighbors = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            # Skip the center tile (the original tilecode)
            if dx == 0 and dy == 0:
                continue
            # Calculate the new x and y
            nx = x + dx
            ny = y + dy
            # Ignore tiles with negative coordinates
            if nx >= 0 and ny >= 0:
                # Add the neighbor's tilecode to the list
                neighbors.append(f"z{z}x{nx}y{ny}")

    return neighbors

def bbox_tilecodes(bbox, zoom):
    """
    Lists all tilecodes intersecting with the bounding box at a specific zoom level.

    Args:
        bbox (list): Bounding box in the format [left, bottom, right, top].
        zoom (int): Zoom level to check.

    Returns:
        list: List of intersecting tilecodes.
    """
    west, south, east, north = bbox
    bbox_geom = box(west, south, east, north)
    
    intersecting_tilecodes = []

    for tile in mercantile.tiles(west, south, east, north, zoom):
        tile_geom = box(*mercantile.bounds(tile))
        if bbox_geom.intersects(tile_geom):
            tilecode = f'z{zoom}x{tile.x}y{tile.y}'
            intersecting_tilecodes.append(tilecode)

    return intersecting_tilecodes

def feature_tilecodes(geometry, zoom):
    """
    Lists all tilecodes intersecting with the Shapely geometry at a specific zoom level.

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The Shapely geometry to check for intersections.
        zoom (int): Zoom level to check.

    Returns:
        list: List of intersecting tilecodes.
    """
    intersecting_tilecodes = []

    for tile in mercantile.tiles(*geometry.bounds, zoom):
        tile_geom = box(*mercantile.bounds(tile))
        if geometry.intersects(tile_geom):
            tilecode = f'z{zoom}x{tile.x}y{tile.y}'
            intersecting_tilecodes.append(tilecode)
    return intersecting_tilecodes
import argparse
import pandas as pd
from shapely.geometry import Polygon, box
from shapely.wkt import loads
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.geocode.geocode2geojson import fix_eaggr_wkt
from pyproj import Geod
import fiona
import math
from tqdm import tqdm  # Import tqdm for progress bars

base_cells = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09',
              '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
eaggr_dggs = Eaggr(Model.ISEA4T)

def haversine(lat1, lon1, lat2, lon2):
    """
    Haversine formula to calculate the distance between two points on the Earth's surface.
    
    Parameters:
        lat1, lon1: Latitude and Longitude of the first point in decimal degrees.
        lat2, lon2: Latitude and Longitude of the second point in decimal degrees.
    
    Returns:
        Distance in meters between the two points.
    """
    R = 6371000  # Radius of the Earth in meters

    # Convert latitude and longitude from degrees to radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(delta_phi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return R * c  # Distance in meters

def isea4t2feature(isea4t):
    """
    Converts a DGGS cell into a GeoJSON-like dictionary with additional metadata.

    Parameters:
        isea4t (str): The DGGS cell ID.

    Returns:
        dict: Feature information including geometry and properties.
    """
    eaggr_cell_shape = DggsShape(DggsCell(isea4t), DggsShapeLocation.ONE_FACE)._shape
    cell_to_shp = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell_shape, ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_eaggr_wkt(cell_to_shp)

    # Convert WKT to Shapely geometry
    cell_polygon = loads(cell_to_shp_fixed)
    
    resolution = len(isea4t) - 2
    # Compute centroid
    cell_centroid = cell_polygon.centroid
    center_lat, center_lon = round(cell_centroid.y, 7), round(cell_centroid.x, 7)
    # Compute area using PyProj Geod
    geod = Geod(ellps="WGS84")
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
    # Compute perimeter using PyProj Geod
    edge_len = abs(geod.geometry_area_perimeter(cell_polygon)[1]) / 3  # Perimeter in meters/3

    edge_len_str = f'{round(edge_len, 2)} m'
    cell_area_str = f'{round(cell_area, 2)} m2'

    if cell_area >= 1000_000:
        edge_len_str = f'{round(edge_len / 1000, 2)} km'
        cell_area_str = f'{round(cell_area / (10**6), 2)} km2'

    return {
        "geometry": cell_polygon,
        "isea4t": isea4t,
        "center_lat": center_lat,
        "center_lon": center_lon,
        "cell_area": cell_area_str,
        "edge_len": edge_len_str,
        "resolution": resolution,
    }

def is_within_bbox(cell_polygon, bbox):
    """
    Check if a cell's polygon intersects with the given bounding box.
    """
    bbox_polygon = box(*bbox)  # Create a Shapely polygon for the bounding box
    return cell_polygon.intersects(bbox_polygon)

def get_bbox_center(bbox):
    """
    Calculate the center of the bounding box.
    """
    min_lon, min_lat, max_lon, max_lat = bbox
    center_lon = (min_lon + max_lon) / 2
    center_lat = (min_lat + max_lat) / 2
    return (center_lat, center_lon)

def get_base_cells_for_bbox(bbox):
    """
    Determine which base cells intersect with the bounding box.
    Returns the base cell with the smallest distance to the center of the bounding box.
    """
    bbox_center = get_bbox_center(bbox)
    closest_base_cell = None
    closest_distance = float('inf')

    for base_cell in base_cells:
        feature = isea4t2feature(base_cell)
        if is_within_bbox(feature["geometry"], bbox):
            # Calculate the distance from the base cell's center to the bounding box center using Haversine formula
            base_cell_center = (feature["center_lat"], feature["center_lon"])
            distance = haversine(bbox_center[0], bbox_center[1], base_cell_center[0], base_cell_center[1])
            if distance < closest_distance:
                closest_distance = distance
                closest_base_cell = base_cell

    return closest_base_cell

def get_cells_for_resolution(base_cell, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution.
    """
    current_cells = [base_cell]
    for _ in range(target_resolution):
        next_cells = []
        for cell in current_cells:
            children = eaggr_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells

def process_cells_within_bbox(cells, bbox, chunk_size=10000):
    """
    Process cells in chunks and filter them by the bounding box.
    """
    total_cells = len(cells)
    rows = []

    # Adding tqdm for progress bar
    for i in tqdm(range(0, total_cells, chunk_size), desc="Processing cells", unit="chunk"):
        chunk = cells[i:i + chunk_size]
        for cell in chunk:
            feature = isea4t2feature(cell)
            if feature and is_within_bbox(feature["geometry"], bbox):
                rows.append(feature)
    return rows

def save_as_shapefile(features, output_file):
    """
    Save features as a Shapefile using pandas and shapely.
    """
    df = pd.DataFrame({
        "geometry": [f["geometry"].wkt for f in features],
        "isea4t": [f["isea4t"] for f in features],
        "center_lat": [f["center_lat"] for f in features],
        "center_lon": [f["center_lon"] for f in features],
        "cell_area": [f["cell_area"] for f in features],
        "edge_len": [f["edge_len"] for f in features],
        "resolution": [f["resolution"] for f in features],
    })

    # Convert geometries to Shapely objects for Fiona compatibility
    df["geometry"] = df["geometry"].apply(loads)

    # Save as Shapefile
    schema = {
        "geometry": "Polygon",
        "properties": {
            "isea4t": "str",
            "center_lat": "float",
            "center_lon": "float",
            "cell_area": "str",
            "edge_len": "str",
            "resolution": "int",
        },
    }

    with fiona.open(output_file, "w", driver="ESRI Shapefile", schema=schema, crs="EPSG:4326") as shp:
        for _, row in df.iterrows():
            shp.write({
                "geometry": row["geometry"].__geo_interface__,
                "properties": {
                    "isea4t": row["isea4t"],
                    "center_lat": row["center_lat"],
                    "center_lon": row["center_lon"],
                    "cell_area": row["cell_area"],
                    "edge_len": row["edge_len"],
                    "resolution": row["resolution"],
                },
            })
    print(f"Shapefile saved as {output_file}.shp")

def main():
    parser = argparse.ArgumentParser(description="Generate DGGS grid within a bounding box and save as Shapefile.")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Desired resolution (e.g., 0, 1, 2, 3).")
    parser.add_argument(
        "-b", "--bbox", type=float, nargs=4, 
        help="Bounding box coordinates (min_lon, min_lat, max_lon, max_lat). If not provided, uses the global bounding box."
    )
    args = parser.parse_args()

    # Set the bounding box
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    # Get the closest base cell for the given bounding box
    closest_base_cell = get_base_cells_for_bbox(bbox)
    if not closest_base_cell:
        print("No base cell found for the bounding box.")
        return

    print(f"Closest base cell: {closest_base_cell}")

    # Generate the cells for the given resolution
    cells = get_cells_for_resolution(closest_base_cell, args.resolution)

    # Process the cells that intersect with the bounding box
    features = process_cells_within_bbox(cells, bbox)

    if features:
        save_as_shapefile(features, f"output_resolution_{args.resolution}")
    else:
        print(f"No features found for resolution {args.resolution} within the bounding box.")

if __name__ == "__main__":
    main()

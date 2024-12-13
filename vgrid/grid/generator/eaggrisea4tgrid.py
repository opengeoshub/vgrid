import argparse
import pandas as pd
from shapely.geometry import Polygon
from shapely.wkt import loads
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.geocode.geocode2geojson import fix_eaggr_wkt
from pyproj import Geod
from tqdm import tqdm  # Import tqdm for progress bars
import fiona 
from shapely import wkt

base_cells = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09',
              '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
eaggr_dggs = Eaggr(Model.ISEA4T)

from shapely.geometry import Polygon
from shapely import wkt

def filter_antimeridian_cells(isea4t_boundary, threshold=-128):
    # Get the coordinates of the polygon
    isea4t_boundary_coordinates = isea4t_boundary.exterior.coords

    # Convert the coordinates to a list of tuples (longitude, latitude)
    lon_lat = [(float(lon), float(lat)) for lon, lat in isea4t_boundary_coordinates]

    # Check if any longitude in the boundary is below the threshold
    if any(lon < threshold for lon, _ in lon_lat):
        # Adjust all longitudes accordingly if any is below the threshold
        adjusted_coords = [(lon - 360 if lon > 0 else lon, lat) for lon, lat in lon_lat]
    else:
        # No adjustment needed, use original coordinates
        adjusted_coords = lon_lat

    # Create a new Polygon with the adjusted coordinates
    adjusted_polygon = Polygon(adjusted_coords)

    return adjusted_polygon
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
    cell_polygon = filter_antimeridian_cells(cell_polygon)
    
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

def get_cells_for_resolution(base_cells, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution.
    """
    current_cells = base_cells
    for _ in range(target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc="Generating child cells", unit="cell"):
            children = eaggr_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells

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
    parser = argparse.ArgumentParser(description="Generate full DGGS grid at a specified resolution.")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Desired resolution (e.g., 0, 1, 2, 3).")
    args = parser.parse_args()

    # Generate the cells for the given resolution
    cells = get_cells_for_resolution(base_cells, args.resolution)

    # Process the cells to create features
    features = []
    for cell in tqdm(cells, desc="Processing cells", unit="cell"):
        features.append(isea4t2feature(cell))

    # Save as Shapefile
    save_as_shapefile(features, f"isea4t_{args.resolution}")

if __name__ == "__main__":
    main()

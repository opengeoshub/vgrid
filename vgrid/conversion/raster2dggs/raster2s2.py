import os, argparse, json
from tqdm import tqdm
import rasterio
from vgrid.utils import s2
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
from vgrid.stats.s2stats import s2_metrics
from vgrid.utils.antimeridian import fix_polygon
from vgrid.generator.settings import geodesic_dggs_metrics, geodesic_dggs_to_feature
from math import cos, radians

def get_nearest_s2_resolution(raster_path):
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        pixel_width = transform.a
        pixel_height = -transform.e
        cell_size = pixel_width*pixel_height
        
        if crs.is_geographic: 
            # Latitude of the raster center
            center_latitude = (src.bounds.top + src.bounds.bottom) / 2
            # Convert degrees to meters
            meter_per_degree_lat = 111_320  # Roughly 1 degree latitude in meters
            meter_per_degree_lon = meter_per_degree_lat * cos(radians(center_latitude))

            pixel_width_m = pixel_width * meter_per_degree_lon
            pixel_height_m = pixel_height * meter_per_degree_lat
            cell_size = pixel_width_m*pixel_height_m    
       
    # Find the nearest s2 resolution by comparing the pixel size to the s2 edge lengths
    nearest_resolution = None
    min_diff = float('inf')
        
    # Check resolutions from 0 to 24
    for res in range(25):
        _, _,avg_area = s2_metrics(res)
        diff = abs(avg_area - cell_size)        
        # If the difference is smaller than the current minimum, update the nearest resolution
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res
    
    return nearest_resolution

def convert_numpy_types(obj):
    """ Recursively convert NumPy types to native Python types """
    if isinstance(obj, np.generic):
        return obj.item()  # Convert numpy types like np.uint8 to native Python int
    elif isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(v) for v in obj]
    else:
        return obj

def raster_to_s2(raster_path, resolution=None):
    # Step 1: Determine the nearest s2 resolution if none is provided
    if resolution is None:
        resolution = get_nearest_s2_resolution(raster_path)
        print(f"Nearest s2 resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    s2_cells = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            lat_lng = s2.LatLng.from_degrees(lat, lon)
            cell_id = s2.CellId.from_lat_lng(lat_lng) # return S2 cell at max level 30
            cell_id = cell_id.parent(resolution) # get S2 cell at resolution
            cell_token = s2.CellId.to_token(cell_id) # get Cell ID Token, shorter than cell_id.id()
            s2_cells.add(cell_token)

    # Sample the raster values at the centroids of the s2 hexagons
    s2_data = []
    
    for s2_token in s2_cells:       
        cell_id = s2.CellId.from_token(s2_token)
        s2_cell = s2.Cell(cell_id)
        vertices = [s2_cell.get_vertex(i) for i in range(4)]
        
        shapely_vertices = []
        for vertex in vertices:
            lat_lng = s2.LatLng.from_point(vertex)  # Convert Point to LatLng
            longitude = lat_lng.lng().degrees  # Access longitude in degrees
            latitude = lat_lng.lat().degrees   # Access latitude in degrees
            shapely_vertices.append((longitude, latitude))

        # Close the polygon by adding the first vertex again
        shapely_vertices.append(shapely_vertices[0])  # Closing the polygon
        # Create a Shapely Polygon
        cell_polygon = fix_polygon(Polygon(shapely_vertices)) # Fix antimeridian
        num_edges = 4
        centroid_lat, centroid_lon,avg_edge_len,cell_area =  geodesic_dggs_metrics(cell_polygon,num_edges)
        
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            s2_data.append({
                "s2": s2_token,
                # "centroid": Point(centroid_lon, centroid_lat),
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    # Create the GeoJSON-like structure
    s2_features = []
    for data in tqdm(s2_data, desc="Resampling", unit=" cells"):
        cell_id = s2.CellId.from_token(data['s2'])
        s2_cell = s2.Cell(cell_id)
        if s2_cell:
            # Get the vertices of the cell (4 vertices for a rectangular cell)
            vertices = [s2_cell.get_vertex(i) for i in range(4)]
            # Prepare vertices in (longitude, latitude) format for Shapely
            shapely_vertices = []
            for vertex in vertices:
                lat_lng = s2.LatLng.from_point(vertex)  # Convert Point to LatLng
                longitude = lat_lng.lng().degrees  # Access longitude in degrees
                latitude = lat_lng.lat().degrees   # Access latitude in degrees
                shapely_vertices.append((longitude, latitude))

            # Close the polygon by adding the first vertex again
            shapely_vertices.append(shapely_vertices[0])  # Closing the polygon
            # Create a Shapely Polygon
            cell_polygon = fix_polygon(Polygon(shapely_vertices)) # Fix antimeridian
            num_edges = 4
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            s2_feature = geodesic_dggs_to_feature("s2",data["s2"],resolution,cell_polygon,num_edges)   
            s2_feature["properties"].update(convert_numpy_types(band_properties) )

            s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
 
       
def main():
    parser = argparse.ArgumentParser(description="Convert Raster to S2 Grid")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default= None, help="Resolution of S2 to be generated [0..24]"
    )


    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    
    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return
    if resolution is not None:
        if resolution < 0 or resolution > 24:
            print(f"Please select a resolution in [0..24] range and try again ")
            return

    s2_geojson = raster_to_s2(raster, resolution)
    geojson_name = os.path.splitext(os.path.basename(raster))[0]
    geojson_path = f"{geojson_name}2s2.geojson"
   
    with open(geojson_path, 'w') as f:
        json.dump(s2_geojson, f)
    
    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()

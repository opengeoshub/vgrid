import os, argparse, json
from tqdm import tqdm
import rasterio
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
import csv
import io
from vgrid.stats.rhealpixstats import rhealpix_metrics
from vgrid.generator.settings import geodesic_dggs_metrics, geodesic_dggs_to_feature
from math import cos, radians
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.conversion.dggs2geojson import rhealpix_cell_to_polygon

def get_nearest_rhealpix_resolution(raster_path):
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
       
    nearest_resolution = None
    min_diff = float('inf')
        
    # Check resolutions from 0 to 15
    for res in range(16):
        _, _,avg_area = rhealpix_metrics(res)
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

def raster2rhealpix(rhealpix_dggs, raster_path, resolution=None, format='geojson'):
    # Step 1: Determine the nearest rhealpix resolution if none is provided
    if resolution is None:
        resolution = get_nearest_rhealpix_resolution(raster_path)
        print(f"Nearest rhealpix resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    rhealpix_ids = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            point = (lon, lat)
            rhealpix_cell = rhealpix_dggs.cell_from_point(resolution, point, plane=False)   
            rhealpix_ids.add(str(rhealpix_cell))

    # Sample the raster values at the centroids of the rhealpix hexagons
    rhealpix_data = []
    
    for rhealpix_id in tqdm(rhealpix_ids, desc="Resampling", unit=" cells"):           
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
        num_edges = 4
        if rhealpix_cell.ellipsoidal_shape() == 'dart':
            num_edges = 3        
        centroid_lat, centroid_lon,avg_edge_len,cell_area =  geodesic_dggs_metrics(cell_polygon,num_edges)
        
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            rhealpix_data.append({
                "rhealpix": rhealpix_id,
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    if format.lower() == 'csv':
        import io
        output = io.StringIO()
        if rhealpix_data:
            writer = csv.DictWriter(output, fieldnames=rhealpix_data[0].keys())
            writer.writeheader()
            writer.writerows(rhealpix_data)
        return output.getvalue()
            
    # Create the GeoJSON-like structure
    rhealpix_features = []
    for data in tqdm(rhealpix_data, desc="Converting to GeoJSON", unit=" cells"):
        rhealpix_uids = (data['rhealpix'][0],) + tuple(map(int, data['rhealpix'][1:]))
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        if rhealpix_cell:
            cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
            num_edges = 4
            if rhealpix_cell.ellipsoidal_shape() == 'dart':
                num_edges = 3
            rhealpix_feature = geodesic_dggs_to_feature("rhealpix",str(rhealpix_cell),resolution,cell_polygon,num_edges)   
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            rhealpix_feature["properties"].update(convert_numpy_types(band_properties))
            rhealpix_features.append(rhealpix_feature)

   
    return {
        "type": "FeatureCollection",
        "features": rhealpix_features
    }

def raster2rhealpix_cli():
    """Command line interface for raster2rhealpix"""
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to rHEALPix DGGS")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default=None, 
        help="Resolution [0..15]"
    )
    parser.add_argument(
        '-f', '--format', type=str, required=False, default='geojson',
        choices=['geojson', 'csv'], help="Output format (geojson or csv)"
    )

    args = parser.parse_args()
    
    E = WGS84_ELLIPSOID
    rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
    
    if not os.path.exists(args.raster):
        raise FileNotFoundError(f"The file {args.raster} does not exist.")
    
    if args.resolution is not None and (args.resolution < 0 or args.resolution > 15):
        raise ValueError("Resolution must be in range [0..15]")

    result = raster2rhealpix(rhealpix_dggs, args.raster, args.resolution, args.format)
    
    # Generate output filename
    base_name = os.path.splitext(os.path.basename(args.raster))[0]
    
    if args.format.lower() == 'csv':
        output_path = f"{base_name}2rhealpix.csv"
        with open(output_path, 'w', newline='') as f:
            f.write(result)
    else:
        output_path = f"{base_name}2rhealpix.geojson"
        with open(output_path, 'w') as f:
            json.dump(result, f)
    
    print(f"Output saved as {output_path}")        

import os, argparse, json
from tqdm import tqdm
import rasterio
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
import csv
from vgrid.stats.isea4tstats import isea4t_metrics
from vgrid.utils.antimeridian import fix_polygon
from vgrid.generator.settings import geodesic_dggs_metrics, geodesic_dggs_to_feature
from math import cos, radians
import platform
if (platform.system() == 'Windows'):   
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.generator.isea4tgrid import fix_isea4t_wkt, fix_isea4t_antimeridian_cells
    from vgrid.generator.settings import isea4t_res_accuracy_dict
from shapely.wkt import loads


def get_nearest_isea4t_resolution(isea4t_dggs,raster_path):
    if (platform.system() == 'Windows'): 
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
            
        # Check resolutions from 0 to 23
        for res in range(24):
            _, _,avg_area,_ = isea4t_metrics(isea4t_dggs,res)
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

def raster2isea4t(isea4t_dggs, raster_path, resolution=None, format='geojson'):
    if (platform.system() == 'Windows'):  
        # Step 1: Determine the nearest isea4t resolution if none is provided
        if resolution is None:
            resolution = get_nearest_isea4t_resolution(isea4t_dggs,raster_path)
            print(f"Nearest isea4t resolution determined: {resolution}")

        # Open the raster file to get metadata and data
        with rasterio.open(raster_path) as src:
            raster_data = src.read()  # Read all bands
            transform = src.transform
            width, height = src.width, src.height
            band_count = src.count  # Number of bands in the raster

        isea4t_ids = set()
        
        for row in range(height):
            for col in range(width):
                lon, lat = transform * (col, row)
                max_accuracy =  isea4t_res_accuracy_dict[39] # maximum cell_id length with 41 characters
                lat_long_point = LatLongPoint(lat, lon, max_accuracy)
                isea4t_cell_max_accuracy = isea4t_dggs.convert_point_to_dggs_cell(lat_long_point)
                cell_id_len = resolution+2
                isea4t_cell = DggsCell(isea4t_cell_max_accuracy._cell_id[:cell_id_len])
                isea4t_ids.add(isea4t_cell._cell_id)

        # Sample the raster values at the centroids of the isea4t hexagons
        isea4t_data = []
        
        for isea4t_id in tqdm(isea4t_ids, desc="Resampling", unit=" cells"):
            cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(isea4t_id),ShapeStringFormat.WKT)
            cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
            if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14')\
                or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
                cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
            
            num_edges = 3
            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))        
            centroid_lat, centroid_lon,avg_edge_len,cell_area =  geodesic_dggs_metrics(cell_polygon,num_edges)            
            col, row = ~transform * (centroid_lon, centroid_lat)
            
            if 0 <= col < width and 0 <= row < height:
                # Get the values for all bands at this centroid
                values = raster_data[:, int(row), int(col)]
                isea4t_data.append({
                    "isea4t": isea4t_id,
                    **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
                })
        
        if format.lower() == 'csv':
            import io
            output = io.StringIO()
            if isea4t_data:
                writer = csv.DictWriter(output, fieldnames=isea4t_data[0].keys())
                writer.writeheader()
                writer.writerows(isea4t_data)
            return output.getvalue()
            
        # Create the GeoJSON-like structure
        isea4t_features = []
        for data in tqdm(isea4t_data, desc="Converting to GeoJSON", unit=" cells"):
            isea4t_id = data['isea4t']
            isea4t_cell = DggsCell(isea4t_id)
            cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(isea4t_cell,ShapeStringFormat.WKT)
            cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
            if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14')\
                or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
                cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
                        
            if cell_to_shape_fixed:
                cell_resolution = len(isea4t_id)-2
                num_edges = 3
                cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))

                isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,cell_resolution,cell_polygon,num_edges)   
                band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
                isea4t_feature["properties"].update(convert_numpy_types(band_properties) )
                isea4t_features.append(isea4t_feature)

        return {
            "type": "FeatureCollection",
            "features": isea4t_features
        }

def raster2isea4t_cli():
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to Open-Eaggr ISEA4T DGGS")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default=None, help="Resolution [0..23]"
    )

    parser.add_argument(
        '-f', '--format', type=str, required=False, default='geojson', 
        choices=['geojson', 'csv'], help="Output format (geojson or csv)"
    )

    if (platform.system() == 'Windows'):
        args = parser.parse_args()
        raster = args.raster
        resolution = args.resolution
        format = args.format
        isea4t_dggs = Eaggr(Model.ISEA4T)

        if not os.path.exists(raster):
            print(f"Error: The file {raster} does not exist.")
            return
        if resolution is not None:
            if resolution < 0 or resolution > 23:
                print(f"Please select a resolution in [0..23] range and try again ")
                return

        result = raster2isea4t(isea4t_dggs, raster, resolution, format)
        base_name = os.path.splitext(os.path.basename(raster))[0]
        
        if format.lower() == 'csv':
            output_path = f"{base_name}2isea4t.csv"
            with open(output_path, 'w', newline='') as f:
                if result:
                    writer = csv.DictWriter(f, fieldnames=result[0].keys())
                    writer.writeheader()
                    writer.writerows(result)
        else:
            output_path = f"{base_name}2isea4t.geojson"
            with open(output_path, 'w') as f:
                json.dump(result, f)
        
        print(f"Output saved as {output_path}")

if __name__ == "__main__":
    raster2isea4t_cli()

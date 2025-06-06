from vgrid.utils import s2, olc, geohash, mercantile
import h3
import argparse
import json
from tqdm import tqdm
import os, re
from vgrid.stats.s2stats import s2_metrics
from vgrid.stats.rhealpixstats import rhealpix_metrics
from vgrid.stats.isea4tstats import isea4t_metrics
from vgrid.stats.qtmstats import qtm_metrics
from vgrid.stats.olcstats import olc_metrics
from vgrid.stats.geohashstats import geohash_metrics
from vgrid.stats.tilecodestats import tilecode_metrics
from vgrid.stats.quadkeystats import quadkey_metrics

from shapely.geometry import shape
from vgrid.generator import h3grid, s2grid, rhealpixgrid, isea4tgrid, qtmgrid, olcgrid, geohashgrid, tilecodegrid, quadkeygrid
from numbers import Number
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from pyproj import Geod
import platform
if (platform.system() == 'Windows'):
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.generator.settings import isea4t_res_accuracy_dict
    

geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID

def get_nearest_resolution(geojson_features, from_dggs, to_dggs, from_field=None):
    if not (from_field):
        from_field = from_dggs    
    try:
        from_dggs_id = geojson_features['features'][0]['properties'][from_field]
    except:
        print(f"There is no valid DGGS IDs found in <{from_field}> field.")
        return 
    try:
        if (from_dggs == 'h3'):
            from_resolution = h3.get_resolution(from_dggs_id)
            from_area = h3.average_hexagon_area(from_resolution,unit='m^2')
    
        elif (from_dggs == 's2'):
            s2_id = s2.CellId.from_token(from_dggs_id)
            from_resolution = s2_id.level()
            _, _, from_area = s2_metrics(from_resolution)        
        
        elif (from_dggs == 'rhealpix'):
            rhealpix_uids = (from_field[0],) + tuple(map(int, from_field[1:]))
            rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3) 
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            from_resolution = rhealpix_cell.resolution   
            _, _, from_area = rhealpix_metrics(from_resolution)        
       
        elif (from_dggs == 'isea4t'):
            if (platform.system() == 'Windows'): 
                from_resolution = len(from_field)-2
                _, _, from_area,_ = isea4t_metrics(isea4t_dggs,from_resolution)        
        
        elif (from_dggs == 'qtm'):
            from_resolution = len(from_field)
            _, _, from_area = qtm_metrics(from_resolution)        
        
        elif (from_dggs == 'olc'):
            coord = olc.decode(from_field)  
            from_resolution = coord.codeLength 
            _, _, from_area = olc_metrics(from_resolution)        
        
        elif (from_dggs == 'geohash'):
            from_resolution = len(from_field)
            _, _, from_area = geohash_metrics(from_resolution)        
        
        elif (from_dggs == 'tilecode'):
            match = re.match(r'z(\d+)x(\d+)y(\d+)', from_field)
            from_resolution = int(match.group(1))
            _, _, from_area = tilecode_metrics(from_resolution)        

        elif (from_dggs == 'quadkey'):
            tile = mercantile.quadkey_to_tile(from_resolution)    
            from_resolution = tile.z
            _, _, from_area = quadkey_metrics(from_resolution)        

    except:
        return         
    
    nearest_resolution = None
    min_diff = float('inf')
   
    if to_dggs == 'h3':                
        for res in range (16):
            avg_area = h3.average_hexagon_area(res,unit='m^2')
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res

    elif to_dggs == 's2':                
        for res in range(31):
            _, _, avg_area = s2_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res
    
    elif to_dggs == 'rhealpix':                
        for res in range(16):
            _, _, avg_area = rhealpix_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res
   
    elif to_dggs == 'isea4t':  
        if (platform.system() == 'Windows'):   
            isea4t_dggs = Eaggr(Model.ISEA4T)   
            for res in range(26):
                _, _, avg_area,_ = isea4t_metrics(isea4t_dggs, res)
                diff = abs(avg_area - from_area)        
                # If the difference is smaller than the current minimum, update the nearest resolution
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res
    
    elif to_dggs == 'qtm':                
        for res in range(1, 25):
            _, _, avg_area = qtm_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res

    elif to_dggs == 'olc':                
        for res in [2,4,6,8,10,11,12,13,14,15]:
            _, _, avg_area = olc_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res
    
    elif to_dggs == 'geohash':                
        for res in range(1, 11):
            _, _, avg_area = geohash_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res
    
    elif to_dggs == 'tilecode':                
        for res in range(30):
            _, _, avg_area = tilecode_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res

    elif to_dggs == 'quadkey':                
        for res in range(30):
            _, _, avg_area = quadkey_metrics(res)
            diff = abs(avg_area - from_area)        
            # If the difference is smaller than the current minimum, update the nearest resolution
            if diff < min_diff:
                min_diff = diff
                nearest_resolution = res

    return nearest_resolution

def generate_grid(geojson_features, to_dggs, resolution):
    # bbox = extract_bbox(geojson_data)  # Extract bounding box from GeoJSON
    dggs_grid = {}
    if to_dggs == 'h3':
        dggs_grid = h3grid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 's2':
        dggs_grid = s2grid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 'rhealpix':
        dggs_grid = rhealpixgrid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 'isea4t':
        if (platform.system() == 'Windows'): 
            isea4t_dggs = Eaggr(Model.ISEA4T)
            dggs_grid = isea4tgrid.generate_grid_resample(isea4t_dggs,resolution, geojson_features)
    elif to_dggs == 'qtm':
        dggs_grid = qtmgrid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 'olc':
        dggs_grid = olcgrid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 'geohash':
        dggs_grid = geohashgrid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 'tilecode':
        dggs_grid = tilecodegrid.generate_grid_resample(resolution, geojson_features)
    elif to_dggs == 'quadkey':
        dggs_grid = quadkeygrid.generate_grid_resample(resolution, geojson_features)

    return dggs_grid

def resampling(layer1, layer2, resample_field):
    try:
        layer1_features = []
        for feature in layer1["features"]:
            if resample_field not in feature["properties"]:
                raise ValueError(f"There is no <{resample_field}> field in the input GeoJSON feattures.")
            geom = shape(feature["geometry"])
            value = feature["properties"][resample_field]
            layer1_features.append((geom, value))
    except ValueError as e:
        print(e)
        return layer2

    resampled_features = []
    
    for feature in tqdm(layer2["features"], desc="Resampling", unit=' cells'):
        layer2_shape = shape(feature["geometry"])
        resampled_value = 0
        intersected_parts = []
        
        for l1_shape, l1_resample in layer1_features:
            if layer2_shape.intersects(l1_shape):
                if not isinstance(l1_resample, Number):
                    print(f"\n Mon-numeric values found in <{resample_field}>. Resampled field calculation failed.")
                    return layer2
               
                intersection = layer2_shape.intersection(l1_shape)
                intersected_parts.append(intersection)
                if not intersection.is_empty:
                    proportion = intersection.area / l1_shape.area
                    resampled_value += l1_resample * proportion
            
            if not intersected_parts:
                continue  # Skip features that have no intersection

        # Add resampled fileld to properties
        feature["properties"][resample_field] = round(resampled_value, 3)
        resampled_features.append(feature)

    return {
        "type": "FeatureCollection",
        "features": resampled_features
    }
        
        

def main():
    parser = argparse.ArgumentParser(description="DGGS Resample")
    dggs_options = ['h3', 's2', 'rhealpix', 'isea4t', 'qtm', 'olc', 'geohash', 'tilecode', 'quadkey']

    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Input DGGS"
    )

    parser.add_argument(
        '-fromdggs', '--fromdggs',
        type=str,
        choices=dggs_options,
        required=True,
        help="Input DGGS. Choose from: " + ", ".join(dggs_options)
    )
    parser.add_argument('-fromfield', '--fromfield', help="Input DGGS ID field")
  
    parser.add_argument(
        '-resamplefield', '--resamplefield', type=str, help="Numeric field for resampling"
    )

    parser.add_argument(
        '-todggs', '--todggs',
        type=str,
        choices=dggs_options,
        required=True,
        help="Output DGGS"
    )

    parser.add_argument(
        '-r', '--resolution', type=int, help="Output resolution"
    )

    args = parser.parse_args()
    geojson = args.geojson
    from_dggs = args.fromdggs
    from_field = args.fromfield
    resample_field = args.resamplefield
    to_dggs = args.todggs
    to_resolution = args.resolution
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    if (from_dggs == to_dggs):
        print("To DGGS must be different with From DGGS")
        return
     
    with open(geojson) as f:
        geojson_features = json.load(f)
    
    if not geojson_features['features']:
        raise ValueError("GeoJSON contains no features.")
    
    if not to_resolution:
        to_resolution = get_nearest_resolution(geojson_features, from_dggs, to_dggs,from_field)
    if to_resolution:
        resampled_features = generate_grid(geojson_features, to_dggs, to_resolution)
        if  resample_field: 
            resampled_features = resampling(geojson_features, resampled_features,resample_field)
            
        # Define the GeoJSON file path
        if resampled_features:
            geojson_path = f"{from_dggs}_to_{to_dggs}_{to_resolution}.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(resampled_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")
        else:
            print('DGGS Resample failed.')
    else:
        print ('There is no appropriate resolutions found in the target DGGS. DGGS Resample failed.')
        
if __name__ == "__main__":
    main()

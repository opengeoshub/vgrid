
# https://github.com/corteva/gars-field
from vgrid.utils import georef
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
import shapely.geometry as geom
from pyproj import CRS

def generate_georef_grid():
    # Define bounds for the whole planet
    lon_min, lon_max = -180.0, 180.0
    lat_min, lat_max = -90.0, 90.0
    
    # Initialize a list to store GARS grid polygons
    gars_grid = []
    res,step = 0,0
    if res == 0:
        step = 15/60
    # Use numpy to generate ranges with floating-point steps
    longitudes = np.arange(lon_min, lon_max, step)
    latitudes = np.arange(lat_min, lat_max, step)
    
    # Loop over longitudes and latitudes in 30-minute intervals
    for lon in longitudes:
        for lat in latitudes:
            # Create a polygon for each GARS cell
            poly = geom.box(lon, lat, lon + step, lat + step)
            georef_code = georef.encode(lat,lon,res) 
            # print(poly)
            gars_grid.append({'geometry': poly, 'georef': str(georef_code)})
    
    # print(gars_grid)
    
    # # Create a GeoDataFrame
    gars_gdf = gpd.GeoDataFrame(gars_grid, crs=CRS.from_epsg(4326))
    
    # # Save the grid
    gars_gdf.to_file(f'./data/georef/georef_grid_{res}.geojson', driver='GeoJSON')

# Run the function
generate_georef_grid()
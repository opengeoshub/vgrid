import geopandas as gpd
from shapely.geometry import Polygon
from pyproj import CRS

def create_utm_grid(minx, miny, maxx, maxy, cell_size, crs):
    # Calculate the number of rows and columns based on cell size
    rows = int((maxy - miny) / cell_size)
    cols = int((maxx - minx) / cell_size)
    
    # Initialize a list to hold grid polygons
    polygons = []
    
    for i in range(cols):
        for j in range(rows):
            # Calculate the bounds of the cell
            x1 = minx + i * cell_size
            x2 = x1 + cell_size
            y1 = miny + j * cell_size
            y2 = y1 + cell_size
            
            # Create the polygon for the cell
            polygons.append(Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)]))
    
    # Create a GeoDataFrame with the polygons and set the CRS
    grid = gpd.GeoDataFrame({'geometry': polygons}, crs=crs)
    
    return grid

# Define the bounding box in UTM coordinates (minx, miny, maxx, maxy) for the Northern Hemisphere
# Example for UTM zone 48N (EPSG:32648)
# bbox = (100000, 0, 900000, 9500000) # for the North 
# bbox = (100000, 100000, 900000, 10000000) # for the South 
bbox = (100000, 0, 900000, 10000000) #  # for both
cell_size = 10000  # Cell size in meters

# Create the grid with UTM CRS
epsg_code = 32648

crs = CRS.from_epsg(epsg_code)
grid = create_utm_grid(*bbox, cell_size, crs)

# Save the grid as a polygon shapefile
grid.to_file(f'utm_grid_{epsg_code}_polygons.shp')

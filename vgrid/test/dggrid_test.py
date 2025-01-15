# import geopandas
import shapely
from shapely import Point
import geopandas as gpd
from vgrid.utils.dggrid4py import DGGRIDv7

# create an inital instance that knows where the dggrid tool lives, configure temp workspace and log/stdout output
dggrid_instance = DGGRIDv7(executable='/usr/local/bin/dggrid', working_dir='.', capture_logs=False, silent=True, tmp_geo_out_legacy=False, debug=False)


# # global ISEA4T grid at resolution 5 into GeoDataFrame to Shapefile
# gdf1 = dggrid_instance.grid_cell_polygons_for_extent('ISEA4T', 4)
# print(gdf1.head())
# gdf1.to_file('isea4t_4.geojson')

# gdf_centroids = dggrid_instance.grid_cell_centroids_for_extent(dggs_type='ISEA7H', resolution=4, mixed_aperture_level=None, clip_geom=None)

# # clip extent
clip_bound = shapely.geometry.box(20.2,57.00, 28.4,60.0 )

dggrid_cell =  dggrid_instance.grid_cell_polygons_from_cellids(['1'],dggs_type='ISEA4T',resolution = 0,split_dateline=True)    
gdf = dggrid_cell.set_geometry("geometry")  # Ensure the geometry column is set

# Check and set CRS to EPSG:4326 if needed
if gdf.crs is None:
    gdf.set_crs(epsg=4326, inplace=True)
elif not gdf.crs.equals("EPSG:4326"):
    gdf = gdf.to_crs(epsg=4326)

# Convert to GeoJSON string
geojson_str = gdf.to_json()

# Save to a GeoJSON file
gdf.to_file("output.geojson", driver="GeoJSON")

print(geojson_str)

# # ISEA7H grid at resolution 9, for extent of provided WGS84 rectangle into GeoDataFrame to Shapefile
# gdf3 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 9, clip_geom=clip_bound)
# print(gdf3.head())
# gdf3.to_file('grids/est_shape_isea7h_9.shp')

# # generate cell and areal statistics for a ISEA7H grids from resolution 0 to 8 (return a pandas DataFrame)
# df1 = dggrid_instance.grid_stats_table('ISEA3H', 8)
# print(df1.head(8))
# df1.to_csv('isea7h_8_stats.csv', index=False)

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return Polygons with cell IDs as GeoDataFrame
# latitude, longitude = 10.775275567242561, 106.70679737574993
# point = Point(longitude, latitude)
# res = 1

# Create a GeoDataFrame
# geodf_points_wgs84 = gpd.GeoDataFrame([{'geometry': point}], crs="EPSG:4326")

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return Polygons with cell IDs as GeoDataFrame
# gdf4 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only = True, dggs_type = 'ISEA7H',resolution = res)
# print(gdf4.head())
# gdf4.to_file('polycells_from_points_isea7h_1.geojson')

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return cell IDs added as column to the points GDF
# gdf5 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only=True, dggs_type='ISEA4H', resolution=8)
# print(gdf5.head())
# gdf5.to_file('geopoint_cellids_from_points_isea4h_8.shp')

# generate DGGS grid cell polygons based on 'cell_id_list' (a list or np.array of provided cell_ids)
# gdf6 = dggrid_instance.grid_cell_polygons_from_cellids(cell_id_list=[1, 4, 8], 'ISEA7H', 5)
# print(gdf6.head())
# gdf6.to_file('from_seqnums_isea7h_5.shp')

# v0.2.6 API update split at dateline for cartesian GIS tools
# gdf7 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 3, split_dateline=True,output_address_type='Z7')
# gdf7.to_file('global_isea7h_3_interrupted.geojson', driver='GeoJSON')

# gdf_z1 = dggrid_instance.grid_cell_polygons_for_extent('IGEO7', 5, clip_geom=clip_bound, output_address_type='Z7_STRING')
# print(gdf_z1.head(3))
# # gdf_z1.to_file('gdf_z1.geojson',driver='GeoJSON')

# df_z1 = dggrid_instance.guess_zstr_resolution(gdf_z1['name'].values, 'IGEO7', input_address_type='Z7_STRING')
# print(df_z1.head(3))

# df_q2di = dggrid_instance.address_transform(gdf_z1['name'].values, 'IGEO7', 5, input_address_type='Z7_STRING', output_address_type='Q2DI')
# print(df_q2di.head(3))

# df_tri = dggrid_instance.address_transform(gdf_z1['name'].values, 'IGEO7', 5, input_address_type='Z7_STRING', output_address_type='PROJTRI')
# print(df_tri.head(3))

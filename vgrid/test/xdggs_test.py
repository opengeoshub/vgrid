# import xarray as xr
# import xdggs

# ds = xdggs.tutorial.open_dataset("air_temperature", "h3")

# # Decode DGGS coordinates
# ds_idx = ds.pipe(xdggs.decode)

# # Assign geographical coordinates
# # ds_idx = ds_idx.dggs.assign_latlon_coords()

# # Interactive visualization
# ds_idx['air'].isel(time=0).compute().dggs.explore(center=0, cmap="viridis", alpha=0.5)

import xarray as xr
import xdggs
import geopandas as gpd
from shapely.geometry import shape, Polygon
_ = xr.set_options(display_expand_data=False)

original_ds = xdggs.tutorial.open_dataset("air_temperature", "healpix").load()
air_temperature = original_ds.drop_vars(["lat", "lon"])

df = original_ds.to_dataframe().reset_index()

# Step 3: Save the DataFrame to a CSV file
output_csv = "air_temperature.csv"
df.to_csv(output_csv, index=False)


ds = air_temperature.pipe(xdggs.decode)
print (ds)


cell_boundaries = ds.dggs.cell_boundaries()
print(cell_boundaries)
polygons = cell_boundaries.values  # Extract the polygon geometries
cell_ids = cell_boundaries.coords["cell_ids"].values  # Extract the cell IDs

# Step 2: Create a GeoDataFrame
# Ensure the polygons are valid Shapely geometries
geometry = [Polygon(p) if isinstance(p, Polygon) else shape(p) for p in polygons]

gdf = gpd.GeoDataFrame(
    {"cell_id": cell_ids},
    geometry=geometry,
    crs="EPSG:4326"  # Assuming WGS84 CRS
)

# Step 3: Export to GeoJSON
output_file = "healpix_cells.geojson"
gdf.to_file(output_file, driver="GeoJSON")

print(f"GeoJSON saved to {output_file}")


# cell_ids = air_temperature["cell_ids"].values  # Extract H3 cell IDs
# print (cell_ids)
# geometries = [Polygon(h3.cell_to_boundary(str(h3_id))) for h3_id in cell_ids]

# # Step 3: Extract the air temperature data and reshape it
# data = air_temperature["air"].mean(dim="time").values  # Average air temperature over time
# df = pd.DataFrame({"cell_id": cell_ids, "air_temperature": data})

# # Step 4: Create a GeoDataFrame
# gdf = gpd.GeoDataFrame(
#     df,
#     geometry=geometries,
#     crs="EPSG:4326"  # Assuming WGS84 CRS
# )

# # Step 5: Export to GeoJSON
# output_file = "air_temperature_cells.geojson"
# gdf.to_file(output_file, driver="GeoJSON")

# print(f"GeoJSON saved to {output_file}")



# print(ds.dggs.sel_latlon([48.0, 48.1], -5.0))


# # Step 2: Generate a DGGS grid covering the whole world
# grid = dggs.grid(bounds=(-180, -90, 180, 90))  # Global bounds

# # Step 3: Convert the DGGS grid to an xarray DataArray
# dggs_xarray = xr.DataArray(
#     data=grid.data,
#     dims=["cell_id"],
#     coords={"cell_id": grid.ids}
# )

# # Display the xarray DataArray structure
# print("DGGS Xarray DataArray:")
# print(dggs_xarray)

# # Step 4: Convert DGGS grid to a GeoDataFrame
# gdf = grid.to_geopandas()

# # Display the GeoDataFrame structure
# print("\nDGGS GeoDataFrame:")
# print(gdf.head())

# # Step 5: Export the GeoDataFrame to GeoJSON
# output_file = "dggs_cells.geojson"
# gdf.to_file(output_file, driver="GeoJSON")

# print(f"\nDGGS cells have been saved to {output_file}")

import geodatasets
import geopandas as gpd
import matplotlib.pyplot as plt
import plotly.express as px
from shapely.geometry import MultiPoint, Point
from shapely.ops import voronoi_diagram
from srai.regionalizers import VoronoiRegionalizer
import pandas as pd

# Function to generate flat Voronoi regions
def generate_flat_voronoi_diagram_regions(
    seeds_gdf: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    points = MultiPoint(seeds_gdf.geometry.values)

    # Generate 2D diagram
    regions = voronoi_diagram(points)

    # Map geometries to GeoDataFrame
    flat_voronoi_regions = gpd.GeoDataFrame(
        geometry=list(regions.geoms),
        crs="EPSG:4326",
    )
    # Apply indexes from the seeds dataframe
    flat_voronoi_regions.index = gpd.pd.Index(
        flat_voronoi_regions.sjoin(seeds_gdf)["index_right"],
        name="region_id",
    )

    # Clip to Earth boundaries
    flat_voronoi_regions.geometry = flat_voronoi_regions.geometry.clip_by_rect(
        xmin=-180, ymin=-90, xmax=180, ymax=90
    )
    return flat_voronoi_regions

# Reading the AED world GeoJSON dataset
aed_world_gdf = gpd.read_file(
    # "https://raw.githubusercontent.com/RaczeQ/medium-articles/main/articles/spherical-geovoronoi/aed_world.geojson"
   './50_cities.geojson'
)
aed_world_gdf_subset = aed_world_gdf.head(50)  # Test with a smaller subset

# Generating flat Voronoi regions for AED world
aed_flat_voronoi_regions = generate_flat_voronoi_diagram_regions(aed_world_gdf)

# Saving the flat Voronoi regions to GeoJSON
aed_flat_voronoi_regions.to_file("flat_voronoi_regions.geojson", driver="GeoJSON")

# Generating spherical Voronoi regions
aed_spherical_voronoi_regions = VoronoiRegionalizer(
    seeds=aed_world_gdf_subset, max_meters_between_points=1_000

).transform()

# batch_size = 10  # Adjust this based on memory capacity
# batches = [aed_world_gdf.iloc[i:i + batch_size] for i in range(0, len(aed_world_gdf), batch_size)]

# all_regions = []
# for batch in batches:
#     regions = VoronoiRegionalizer(
#         seeds=batch, 
#         max_meters_between_points=1_000
#     ).transform()
#     all_regions.append(regions)

# # Combine all batch results into a single GeoDataFrame
# full_voronoi_regions = gpd.GeoDataFrame(pd.concat(all_regions, ignore_index=True))


# # Saving the spherical Voronoi regions to GeoJSON
aed_spherical_voronoi_regions.to_file("spherical_voronoi_regions.geojson", driver="GeoJSON")

# You can visualize them if you like
aed_spherical_voronoi_regions.plot()
plt.show()

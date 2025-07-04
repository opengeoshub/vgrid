import shapely
from pyproj import Geod
from shapely.geometry import Polygon, Point, LineString
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

# Initialize Geod with WGS84 ellipsoid
geod = Geod(ellps="WGS84")

def shortest_point_distance(points):
    """
    Calculate distances between points in a Shapely geometry.
    If there's only one point, return 0.
    If there are multiple points, calculate Delaunay triangulation and return distances.

    Args:
        points: Shapely Point or MultiPoint geometry

    Returns:
        tuple: shortest_distance
    """
    # Handle single Point
    if isinstance(points, Point):
        return 0  # Single point has no distance to other points

    # Handle MultiPoint with single point
    if len(points.geoms) == 1:
        return 0

    # Generate Delaunay triangulation
    delaunay = shapely.delaunay_triangles(points, only_edges=True)

    # Find the shortest edge
    shortest_distance = float("inf")

    for line in delaunay.geoms:
        # Get the coordinates of the line endpoints
        coords = list(line.coords)
        lon1, lat1 = coords[0]
        lon2, lat2 = coords[1]

        # Calculate the distance in meters using pyproj Geod
        distance = geod.inv(lon1, lat1, lon2, lat2)[2]  # [2] gives the distance in meters
        if distance < shortest_distance:
            shortest_distance = distance

    return shortest_distance if shortest_distance != float("inf") else 0


def shortest_polyline_distance(polylines):
    """
    Calculate the shortest distance between polylines using GeoPandas shortest_line() method.
    If there's only one polyline, return 0.
    If there are multiple polylines, use shortest_line() and return the shortest distance in meters.

    Args:
        polylines: Shapely LineString or MultiLineString geometry, or GeoSeries of LineStrings
        
    Returns:
        float: shortest_distance between polylines in meters
    """
    # Handle single LineString
    if isinstance(polylines, LineString):
        return 0  # Single polyline has no distance to other polylines
    
    # Handle MultiLineString with single line
    if hasattr(polylines, 'geoms') and len(polylines.geoms) == 1:
        return 0
    
    # Handle GeoSeries
    if hasattr(polylines, 'iloc'):
        # Already a GeoSeries
        line_list = list(polylines.geometry)
        gs = polylines
    else:
        # Handle MultiLineString or list
        line_list = list(polylines.geoms) if hasattr(polylines, 'geoms') else [polylines]
        gs = gpd.GeoSeries(line_list)
    
    if len(line_list) < 2:
        return 0
    
    # Calculate shortest distance between all pairs of polylines using shortest_line()
    shortest_distance = float("inf")
    
    for i in range(len(line_list)):
        for j in range(i+1, len(line_list)):
            line1 = line_list[i]
            line2 = line_list[j]
            
            # Check if polylines are disjoint
            if line1.disjoint(line2):
                try:
                    # Create GeoSeries for shortest_line calculation
                    gs1 = gpd.GeoSeries([line1])
                    gs2 = gpd.GeoSeries([line2])
                    
                    # Get shortest line using GeoPandas method
                    shortest_line = gs1.shortest_line(gs2, align=False).iloc[0]
                    
                    if shortest_line and shortest_line.length > 0:
                            # Get the endpoints of the shortest line
                            coords = list(shortest_line.coords)
                            if len(coords) >= 2:
                                lon1, lat1 = coords[0]
                                lon2, lat2 = coords[1]
                                
                                # Calculate geodesic distance in meters
                                distance = geod.inv(lon1, lat1, lon2, lat2)[2]  # [2] gives distance in meters
                                
                                if distance < shortest_distance:
                                    shortest_distance = distance
                except Exception as e:
                    print(f"Error calculating distance between polylines {i} and {j}: {e}")
                    continue
    
    return shortest_distance if shortest_distance != float("inf") else 0 


def shortest_polygon_distance(polygons):
    """
    Calculate the shortest distance between polygons using GeoPandas shortest_line() method.
    If there's only one polygon, return 0.
    If there are multiple polygons, use shortest_line() and return the shortest distance in meters.

    Args:
        polygons: Shapely Polygon or MultiPolygon geometry, or GeoSeries of Polygons
        
    Returns:
        float: shortest_distance between polygons in meters
    """
    # Handle single Polygon
    if isinstance(polygons, Polygon):
        return 0  # Single polygon has no distance to other polygons
    
    # Handle MultiPolygon with single polygon
    if hasattr(polygons, 'geoms') and len(polygons.geoms) == 1:
        return 0
    
    # Handle GeoSeries
    if hasattr(polygons, 'iloc'):
        # Already a GeoSeries
        polygon_list = list(polygons.geometry)
        gs = polygons
    else:
        # Handle MultiPolygon or list
        polygon_list = list(polygons.geoms) if hasattr(polygons, 'geoms') else [polygons]
        gs = gpd.GeoSeries(polygon_list)
    
    if len(polygon_list) < 2:
        return 0
    
    # Calculate shortest distance between all pairs of polygons using shortest_line()
    shortest_distance = float("inf")
    
    for i in range(len(polygon_list)):
        for j in range(i+1, len(polygon_list)):
            polygon1 = polygon_list[i]
            polygon2 = polygon_list[j]
            
            # Check if polygons are disjoint
            if polygon1.disjoint(polygon2):
                try:
                    # Create GeoSeries for shortest_line calculation
                    gs1 = gpd.GeoSeries([polygon1])
                    gs2 = gpd.GeoSeries([polygon2])
                    
                    # Get shortest line using GeoPandas method
                    shortest_line = gs1.shortest_line(gs2, align=False).iloc[0]
                    
                    if shortest_line and shortest_line.length > 0:
                            # Get the endpoints of the shortest line
                            coords = list(shortest_line.coords)
                            if len(coords) >= 2:
                                lon1, lat1 = coords[0]
                                lon2, lat2 = coords[1]
                                
                                # Calculate geodesic distance in meters
                                distance = geod.inv(lon1, lat1, lon2, lat2)[2]  # [2] gives distance in meters
                                
                                if distance < shortest_distance:
                                    shortest_distance = distance
                except Exception as e:
                    print(f"Error calculating distance between polygons {i} and {j}: {e}")
                    continue
    
    return shortest_distance if shortest_distance != float("inf") else 0


def geodesic_distance(
    lat: float, lon: float, length_meter: float
) -> tuple[float, float]:
    """
    Convert meters to approximate degree offsets at a given location.

    Parameters:
        lat (float): Latitude of the reference point
        lon (float): Longitude of the reference point
        length_meter (float): Distance in meters

    Returns:
        (delta_lat_deg, delta_lon_deg): Tuple of degree offsets in latitude and longitude
    """
    # Move north for latitude delta
    lon_north, lat_north, _ = geod.fwd(lon, lat, 0, length_meter)
    delta_lat = lat_north - lat

    # Move east for longitude delta
    lon_east, lat_east, _ = geod.fwd(lon, lat, 90, length_meter)
    delta_lon = lon_east - lon

    return delta_lat, delta_lon


def geodesic_buffer(polygon, distance):
    """
    Create a geodesic buffer around a polygon using pyproj Geod.

    Args:
        polygon: Shapely Polygon geometry
        distance: Buffer distance in meters

    Returns:
        Shapely Polygon: Buffered polygon
    """
    buffered_coords = []
    for lon, lat in polygon.exterior.coords:
        # Generate points around the current vertex to approximate a circle
        circle_coords = [
            geod.fwd(lon, lat, azimuth, distance)[
                :2
            ]  # Forward calculation: returns (lon, lat, back_azimuth)
            for azimuth in range(0, 360, 10)  # Generate points every 10 degrees
        ]
        buffered_coords.append(circle_coords)

    # Flatten the list of buffered points and form a Polygon
    all_coords = [coord for circle in buffered_coords for coord in circle]
    return Polygon(all_coords).convex_hull


def check_predicate(cell_polygon, input_geometry, predicate=None):
    """
    Determine whether to keep an H3 cell based on its relationship with the input geometry.

    Args:
        cell_polygon: Shapely Polygon representing the H3 cell
        input_geometry: Shapely geometry (Polygon, LineString, etc.)
        predicate (str or int): Spatial predicate to apply:
            String values:
                None or "intersects": intersects (default)
                "within": within
                "centroid_within": centroid_within
                "largest_overlap": intersection >= 50% of cell area
            Integer values (for backward compatibility):
                None or 0: intersects (default)
                1: within
                2: centroid_within
                3: intersection >= 50% of cell area

    Returns:
        bool: True if cell should be kept, False otherwise
    """
    # Handle string predicates
    if isinstance(predicate, str):
        predicate_lower = predicate.lower()
        if predicate_lower in ["intersects", "intersect"]:
            return cell_polygon.intersects(input_geometry)
        elif predicate_lower == "within":
            return cell_polygon.within(input_geometry)
        elif predicate_lower in ["centroid_within", "centroid"]:
            return cell_polygon.centroid.within(input_geometry)
        elif predicate_lower in ["largest_overlap", "overlap", "majority"]:
            # intersection >= 50% of cell area
            if cell_polygon.intersects(input_geometry):
                intersection_geom = cell_polygon.intersection(input_geometry)
                if intersection_geom and intersection_geom.area > 0:
                    intersection_area = intersection_geom.area
                    cell_area = cell_polygon.area
                    return (intersection_area / cell_area) >= 0.5
            return False
        else:
            # Unknown string predicate, default to intersects
            return cell_polygon.intersects(input_geometry)

    # Handle integer predicates (backward compatibility)
    elif isinstance(predicate, int):
        if predicate == 0:
            # Default: intersects
            return cell_polygon.intersects(input_geometry)
        elif predicate == 1:
            # within
            return cell_polygon.within(input_geometry)
        elif predicate == 2:
            # centroid_within
            return cell_polygon.centroid.within(input_geometry)
        elif predicate == 3:
            # intersection >= 50% of cell area
            if cell_polygon.intersects(input_geometry):
                intersection_geom = cell_polygon.intersection(input_geometry)
                if intersection_geom and intersection_geom.area > 0:
                    intersection_area = intersection_geom.area
                    cell_area = cell_polygon.area
                    return (intersection_area / cell_area) >= 0.5
            return False
        else:
            # Unknown predicate, default to intersects
            return cell_polygon.intersects(input_geometry)

    else:
        # None or other types, default to intersects
        return cell_polygon.intersects(input_geometry)


#### Test - Shortest Lines Between Polygons
# Example demonstrating shortest lines between polygons

# Create test polygons (3 separated polygons)
# polygons = [
#     Polygon([(0, 0), (2, 2), (0, 2)]),
#     Polygon([(4, 4), (7, 8), (8, 0)]),
#     Polygon([(8, 8), (10, 12), (11, 10)]),
# ]

# # Create GeoSeries
# s = gpd.GeoSeries(polygons,)

# # Calculate shortest lines between all pairs using pairwise approach
# print("Calculating shortest lines between polygons...")
# shortest_lines = s.shortest_line(s, align=False)
# print(shortest_lines)
# # # Create the plot
# fig, ax = plt.subplots(figsize=(10, 8))

# # Plot polygons with different colors
# colors = ['lightblue', 'lightgreen', 'lightcoral']
# for i, (polygon, color) in enumerate(zip(polygons, colors)):
#     x, y = polygon.exterior.xy
#     ax.fill(x, y, alpha=0.7, fc=color, ec='black', linewidth=2, label=f'Polygon {i+1}')
    
#     # Add polygon label
#     centroid = polygon.centroid
#     ax.annotate(f'P{i+1}', (centroid.x, centroid.y), 
#                fontsize=12, ha='center', va='center', weight='bold')

# # Plot shortest lines between all pairs
# line_colors = ['red', 'blue', 'green']
# line_idx = 0

# for i in range(len(polygons)):
#     for j in range(len(polygons)):
#         if i != j:  # Don't connect polygon to itself
#             try:
#                 # Create GeoSeries with proper indices to avoid alignment warning
#                 source_series = gpd.GeoSeries([s.iloc[i]], index=[0])
#                 target_series = gpd.GeoSeries([s.iloc[j]], index=[0])
                
#                 # Get the shortest line with explicit alignment
#                 shortest_lines = source_series.shortest_line(target_series, align=True)
#                 line = shortest_lines.iloc[0]
                
#                 if line and not line.is_empty:
#                     x, y = line.xy
#                     color = line_colors[line_idx % len(line_colors)]
#                     ax.plot(x, y, color=color, linewidth=2, linestyle='--', 
#                            label=f'P{i+1} to P{j+1}' if line_idx < 3 else "")
                    
#                     # Add distance label at midpoint
#                     midpoint = line.interpolate(0.5, normalized=True)
#                     distance = line.length
#                     ax.annotate(f'{distance:.2f}', (midpoint.x, midpoint.y), 
#                                fontsize=8, ha='center', va='bottom',
#                                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
#                     line_idx += 1
                    
#                     print(f"P{i+1} to P{j+1}: Length = {distance:.4f}")
#             except Exception as e:
#                 print(f"Error calculating line from P{i+1} to P{j+1}: {e}")

# # Customize the plot
# ax.set_xlim(-1, 12)
# ax.set_ylim(-1, 13)
# ax.set_xlabel('X Coordinate')
# ax.set_ylabel('Y Coordinate')
# ax.set_title('3 Separated Polygons with Shortest Connecting Lines')
# ax.grid(True, alpha=0.3)
# ax.legend()
# ax.set_aspect('equal')

# plt.tight_layout()
# plt.show()

# print("\nShortest lines calculation completed!")

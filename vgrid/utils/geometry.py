import shapely
from pyproj import Geod
from shapely.geometry import Polygon, Point, LineString
import numpy as np

# Initialize Geod with WGS84 ellipsoid
geod = Geod(ellps="WGS84")


def calculate_point_distances(points):
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
        distance = geod.inv(lon1, lat1, lon2, lat2)[
            2
        ]  # [2] gives the distance in meters
        if distance < shortest_distance:
            shortest_distance = distance

    return shortest_distance


geod = Geod(ellps="WGS84")


def densify_line(line, segment_length):
    total_length = line.length
    if total_length == 0:
        # Degenerate line, just return start and end twice
        coords = list(line.coords)
        if len(coords) == 1:
            coords = coords * 2
        return LineString(coords)
    num_segments = max(1, int(np.ceil(total_length / segment_length)))
    distances = np.linspace(0, total_length, num_segments + 1)
    points = [line.interpolate(d) for d in distances]
    # Ensure at least two points
    if len(points) < 2:
        points = [line.interpolate(0), line.interpolate(total_length)]
    return LineString(points)


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


#### Test
# Read points from GeoJSON file
# with open('./data/shape/point.geojson', 'r', encoding='utf-8') as f:
#     geojson_data = json.load(f)

# # Extract coordinates from the GeoJSON features
# coordinates = []
# for feature in geojson_data['features']:
#     coords = feature['geometry']['coordinates']
#     coordinates.append((coords[0], coords[1]))  # (longitude, latitude)

# print(f"Number of points: {len(coordinates)}")

# # Create MultiPoint from the coordinates
# points = MultiPoint(coordinates)

# # Calculate distances using the function
# shortest_distance = calculate_point_distances(points)
# print(f"Shortest distance: {shortest_distance:.2f} meters")

# # Plot the points
# x_coords = [point.x for point in points.geoms]
# y_coords = [point.y for point in points.geoms]
# plt.scatter(x_coords, y_coords, color='red', s=100, zorder=3)

# # Generate Delaunay triangulation for plotting
# if len(points.geoms) > 1:
#     delaunay = shapely.delaunay_triangles(points, only_edges=True)

#     # Find the shortest edge for highlighting
#     shortest_edge = None
#     for line in delaunay.geoms:
#         coords = list(line.coords)
#         lon1, lat1 = coords[0]
#         lon2, lat2 = coords[1]
#         distance = geod.inv(lon1, lat1, lon2, lat2)[2]
#         if abs(distance - shortest_distance) < 0.01:  # Small tolerance for floating point comparison
#             shortest_edge = line
#             break

#     # Plot the delaunay triangulation edges
#     for line in delaunay.geoms:
#         x_coords = [coord[0] for coord in line.coords]
#         y_coords = [coord[1] for coord in line.coords]

#         # Check if this is the shortest edge
#         if line == shortest_edge:
#             plt.plot(x_coords, y_coords, 'g-', linewidth=3, label='Shortest Edge')
#         else:
#             plt.plot(x_coords, y_coords, 'b-', linewidth=1)

# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.title(f'Delaunay Triangulation')
# plt.grid(True, alpha=0.3)
# plt.axis('equal')
# if len(points.geoms) > 1:
#     plt.legend()
# plt.show()

import re
import json
import os
import requests
from urllib.parse import urlparse
from vgrid.utils.io import (
    validate_dggal_type,
    validate_dggrid_type,
    init_dggs_band_accumulator,
    update_dggs_band_accumulator,
)
import shapely.geometry
import shapely
from pyproj import CRS as PyprojCRS, Geod
from shapely.geometry import Point, Polygon, mapping, MultiPolygon, box
from shapely.wkt import loads
from vgrid.dggs.rhealpixdggs.utils import my_round
from vgrid.dggs.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.dggs.eaggr.eaggr import Eaggr
from vgrid.dggs.eaggr.enums.model import Model
from vgrid.utils.constants import (
    DGGAL_TYPES,
    AUTHALIC_RADIUS,
    ICOSA_EDGE_M,
    STANDARD_METERS_PER_PIXEL,
)
import platform
import math
from dggal import *
import pandas as pd
import numpy as np
from shapely.ops import transform
import pyproj
import geopandas as gpd
from vgrid.utils.constants import AUTHALIC_AREA, DGGS_TYPES, GEOREF_RESOLUTION_DEGREES
from vgrid.utils.constants import DGGS_TYPES
from vgrid.dggs import maidenhead

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.enums.model import Model

    isea4t_dggs = Eaggr(Model.ISEA4T)
    isea3h_dggs = Eaggr(Model.ISEA3H)


from dggal import *

app = Application(appGlobals=globals())
pydggal_setup(app)

# Initialize Geod with WGS84 ellipsoid
geod = Geod(ellps="WGS84")
# geod = Geod(a=6371007.181, f=0)  # sphere


def rhealpix_cell_to_polygon(cell):
    vertices = [
        tuple(my_round(coord, 14) for coord in vertex)
        for vertex in cell.vertices(plane=False)
    ]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    return Polygon(vertices)


def fix_isea4t_wkt(isea4t_wkt):
    coords_section = isea4t_wkt[isea4t_wkt.index("((") + 2 : isea4t_wkt.index("))")]
    coords = coords_section.split(",")
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    fixed_coords = ", ".join(coords)
    return f"POLYGON (({fixed_coords}))"


def check_crossing(lon1: float, lon2: float, validate: bool = True):
    """
    Assuming a minimum travel distance between two provided longitude coordinates,
    checks if the 180th meridian (antimeridian) is crossed.
    """
    if validate and any(abs(x) > 180.0 for x in [lon1, lon2]):
        raise ValueError("longitudes must be in degrees [-180.0, 180.0]")
    return abs(lon2 - lon1) > 180.0


def check_crossing_geom(geom):
    """
    Check if a geometry crosses the antimeridian (180th meridian).

    Args:
        geom: Shapely geometry (Polygon, MultiPolygon)

    Returns:
        bool: True if any part of the geometry crosses the antimeridian, False otherwise
    """
    crossed = False
    try:
        # Handle None or empty geometries
        if geom is None or geom.is_empty:
            return False

        # Handle multi-geometries
        if hasattr(geom, "geoms"):
            # MultiPolygon
            for sub_geom in geom.geoms:
                if check_crossing_geom(sub_geom):
                    crossed = True
                    break
            return crossed

        # Handle single Polygon
        if geom.geom_type == "Polygon":
            # Check if exterior has coordinates
            if geom.exterior is None or len(geom.exterior.coords) == 0:
                return False

            # Check exterior ring only
            p_init = geom.exterior.coords[0]
            for p in range(1, len(geom.exterior.coords)):
                px = geom.exterior.coords[p]
                try:
                    if check_crossing(p_init[0], px[0]):
                        crossed = True
                        break
                except ValueError:
                    crossed = True
                    break
                p_init = px
    except Exception:
        return False

    return crossed


# -----------------------------------------------------------
# Shift to WEST side (only when crossing)
# -----------------------------------------------------------
def shift_west(geom, threshold=-130):
    if not check_crossing_geom(geom):
        return geom  # do nothing

    if hasattr(geom, "geoms"):
        # MultiPolygon - process each sub-geometry
        shifted_geoms = [shift_west(sub_geom, threshold) for sub_geom in geom.geoms]
        return MultiPolygon(shifted_geoms)

    # Single Polygon - extract exterior coordinates
    lon_lat = [(float(lon), float(lat)) for lon, lat in geom.exterior.coords]
    # Only apply shift if any longitude is below threshold
    if any(lon < threshold for lon, _ in lon_lat):
        adjusted_coords = [(lon - 360 if lon > 0 else lon, lat) for lon, lat in lon_lat]
    else:
        adjusted_coords = lon_lat

    # Reconstruct polygon and ensure validity
    shifted_polygon = Polygon(adjusted_coords)
    # if not shifted_polygon.is_valid:
    #     shifted_polygon = make_valid(shifted_polygon)

    return shifted_polygon


# -----------------------------------------------------------
# Shift to EAST side (only when crossing)
# -----------------------------------------------------------
def shift_east(geom, threshold=146):
    if not check_crossing_geom(geom):
        return geom  # do nothing

    if hasattr(geom, "geoms"):
        # MultiPolygon - process each sub-geometry
        shifted_geoms = [shift_east(sub_geom, threshold) for sub_geom in geom.geoms]
        return MultiPolygon(shifted_geoms)

    # Single Polygon - extract exterior coordinates
    lon_lat = [(float(lon), float(lat)) for lon, lat in geom.exterior.coords]
    # Only apply shift if any longitude is above the absolute threshold (mirrored)
    if any(lon > threshold for lon, _ in lon_lat):
        adjusted_coords = [(lon + 360 if lon < 0 else lon, lat) for lon, lat in lon_lat]
    else:
        adjusted_coords = lon_lat

    # Reconstruct polygon and ensure validity
    shifted_polygon = Polygon(adjusted_coords)
    # if not shifted_polygon.is_valid:
    #     shifted_polygon = make_valid(shifted_polygon)

    return shifted_polygon


# -----------------------------------------------------------
# Balanced shift based on centroid location
# -----------------------------------------------------------
def shift_balanced(geom, threshold_west=-130, threshold_east=146):
    """
    Shift geometry based on centroid location:
    - If centroid longitude < -180: apply west shift
    - If centroid longitude >= 180: apply east shift
    - Otherwise: determine based on centroid position

    Args:
        geom: Shapely geometry (Polygon, MultiPolygon)
        threshold: float, optional
            Longitudinal threshold (in degrees) used to detect whether adjustment is needed.
            Defaults to -128.

    Returns:
        Shapely geometry: shifted geometry
    """
    if not check_crossing_geom(geom):
        return geom  # do nothing if not crossing

    if hasattr(geom, "geoms"):
        # MultiPolygon - process each sub-geometry
        shifted_geoms = [
            shift_balanced(sub_geom, threshold_west, threshold_east)
            for sub_geom in geom.geoms
        ]
        return MultiPolygon(shifted_geoms)

    # Get centroid longitude
    centroid = geom.centroid
    centroid_lon = centroid.x

    # Apply shift based on centroid longitude
    if centroid_lon < -180:
        return shift_west(geom, threshold_west)
    elif centroid_lon >= 180:
        return shift_east(geom, threshold_east)
    else:
        # For geometries with centroid in [-180, 180], shift based on sign
        # Negative (western hemisphere) -> west shift, positive -> east shift
        if centroid_lon < 0:
            return shift_west(geom, threshold_west)
        else:
            return shift_east(geom, threshold_east)


def isea4t_cell_to_polygon(isea4t_cell):
    cell_to_shp = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(
        isea4t_cell, ShapeStringFormat.WKT
    )
    cell_to_shp_fixed = fix_isea4t_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)
    return cell_polygon


def get_ease_resolution(ease_id):
    """Get the resolution level of an EASE cell ID."""
    try:
        match = re.match(r"L(\d+)\.(.+)", ease_id)
        if not match:
            raise ValueError(f"Invalid EASE ID format: {ease_id}")
        return int(match.group(1))
    except Exception as e:
        raise ValueError(f"Invalid EASE ID <{ease_id}> : {e}")


def isea3h_cell_to_polygon(isea3h_cell):
    cell_to_shape = isea3h_dggs.convert_dggs_cell_outline_to_shape_string(
        isea3h_cell, ShapeStringFormat.WKT
    )
    cell_to_shp_fixed = fix_isea4t_wkt(cell_to_shape)
    cell_polygon = loads(cell_to_shp_fixed)
    return cell_polygon


def fix_eaggr_wkt(eaggr_wkt):
    coords_section = eaggr_wkt[eaggr_wkt.index("((") + 2 : eaggr_wkt.index("))")]
    coords = coords_section.split(",")
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    fixed_coords = ", ".join(coords)
    return f"POLYGON (({fixed_coords}))"


def graticule_dggs_metrics(cell_polygon):
    min_lon, min_lat, max_lon, max_lat = cell_polygon.bounds
    center_lat = (min_lat + max_lat) / 2
    center_lon = (min_lon + max_lon) / 2
    cell_width = geod.line_length([min_lon, max_lon], [min_lat, min_lat])
    cell_height = geod.line_length([min_lon, min_lon], [min_lat, max_lat])
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
    cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
    return center_lat, center_lon, cell_width, cell_height, cell_area, cell_perimeter


def geodesic_dggs_metrics(cell_polygon, num_edges):
    cell_centroid = cell_polygon.centroid
    center_lat = cell_centroid.y
    center_lon = cell_centroid.x
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
    cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
    avg_edge_len = cell_perimeter / num_edges
    return center_lat, center_lon, avg_edge_len, cell_area, cell_perimeter


def graticule_dggs_to_feature(dggs_type, cell_id, resolution, cell_polygon):
    center_lat, center_lon, cell_width, cell_height, cell_area, cell_perimeter = (
        graticule_dggs_metrics(cell_polygon)
    )
    feature = {
        "type": "Feature",
        "geometry": mapping(cell_polygon),
        "properties": {
            f"{dggs_type}": str(cell_id),
            "resolution": resolution,
            "center_lat": center_lat,
            "center_lon": center_lon,
            "cell_width": cell_width,
            "cell_height": cell_height,
            "cell_area": cell_area,
            "cell_perimeter": cell_perimeter,
        },
    }
    return feature


def geodesic_dggs_to_feature(dggs_type, cell_id, resolution, cell_polygon, num_edges):
    center_lat, center_lon, avg_edge_len, cell_area, cell_perimeter = (
        geodesic_dggs_metrics(cell_polygon, num_edges)
    )
    feature = {
        "type": "Feature",
        "geometry": mapping(cell_polygon),
        "properties": {
            f"{dggs_type}": str(cell_id),
            "resolution": resolution,
            "center_lat": center_lat,
            "center_lon": center_lon,
            "avg_edge_len": avg_edge_len,
            "cell_area": cell_area,
            "cell_perimeter": cell_perimeter,
        },
    }
    return feature


def graticule_dggs_to_geoseries(dggs_type, cell_id, resolution, cell_polygon):
    center_lat, center_lon, cell_width, cell_height, cell_area, cell_perimeter = (
        graticule_dggs_metrics(cell_polygon)
    )
    return {
        f"{dggs_type}": str(cell_id),
        "resolution": resolution,
        "center_lat": center_lat,
        "center_lon": center_lon,
        "cell_width": cell_width,
        "cell_height": cell_height,
        "cell_area": cell_area,
        "cell_perimeter": cell_perimeter,
        "geometry": cell_polygon,
    }


def geodesic_dggs_to_geoseries(dggs_type, cell_id, resolution, cell_polygon, num_edges):
    center_lat, center_lon, avg_edge_len, cell_area, cell_perimeter = (
        geodesic_dggs_metrics(cell_polygon, num_edges)
    )
    return {
        f"{dggs_type}": str(cell_id),
        "resolution": resolution,
        "center_lat": center_lat,
        "center_lon": center_lon,
        "avg_edge_len": avg_edge_len,
        "cell_area": cell_area,
        "cell_perimeter": cell_perimeter,
        "geometry": cell_polygon,
    }


def raster_bbox_wgs84(raster_path: str) -> list[float]:
    """Return raster bounds as ``[minx, miny, maxx, maxy]`` in EPSG:4326."""
    import rasterio
    from pyproj import Transformer

    with rasterio.open(raster_path) as src:
        if src.crs is None:
            raise ValueError(
                "Raster CRS is undefined. Raster conversion requires a valid CRS."
            )
        b = src.bounds
        if src.crs.to_epsg() == 4326:
            return [b.left, b.bottom, b.right, b.top]
        transformer = Transformer.from_crs(src.crs, "EPSG:4326", always_xy=True)
        xs, ys = transformer.transform(
            [b.left, b.left, b.right, b.right],
            [b.bottom, b.top, b.bottom, b.top],
        )
        return [min(xs), min(ys), max(xs), max(ys)]


def read_pixel_centroids(raster_path: str):
    """
    Read valid raster pixels as point geometries (cell centers) in EPSG:4326.

    Returns
    -------
    tuple
        ``(GeoDataFrame, band_count)`` with columns ``band_1``, ``band_2``, …
    """
    import rasterio
    import geopandas as gpd
    from pyproj import Transformer
    from tqdm import tqdm

    records = []
    band_count = 0

    with rasterio.open(raster_path) as src:
        if src.crs is None:
            raise ValueError(
                "Raster CRS is undefined. Raster conversion requires a valid CRS."
            )
        transform = src.transform
        band_count = src.count
        transformer = None
        if src.crs.to_epsg() != 4326:
            transformer = Transformer.from_crs(src.crs, "EPSG:4326", always_xy=True)

        block_height, block_width = src.block_shapes[0]
        total_blocks = ((src.height + block_height - 1) // block_height) * (
            (src.width + block_width - 1) // block_width
        )
        for _, window in tqdm(
            src.block_windows(1),
            total=total_blocks,
            desc="Reading raster pixel centers",
            unit=" block",
        ):
            block = src.read(window=window, masked=True)
            if block.size == 0:
                continue

            block_mask = np.ma.getmaskarray(block)
            valid_pixels = ~np.all(block_mask, axis=0)
            if not np.any(valid_pixels):
                continue

            rows_local, cols_local = np.where(valid_pixels)
            rows_global = rows_local + window.row_off
            cols_global = cols_local + window.col_off

            xs, ys = rasterio.transform.xy(
                transform, rows_global, cols_global, offset="center"
            )
            if transformer is not None:
                lons, lats = transformer.transform(xs, ys)
            else:
                lons, lats = xs, ys

            for idx, (lon, lat) in enumerate(zip(lons, lats)):
                pixel_vals = block[:, rows_local[idx], cols_local[idx]]
                rec = {"geometry": shapely.geometry.Point(float(lon), float(lat))}
                for bi in range(band_count):
                    v = pixel_vals[bi]
                    rec[f"band_{bi + 1}"] = None if np.ma.is_masked(v) else float(v)
                records.append(rec)

    if not records:
        return gpd.GeoDataFrame(columns=["geometry"], crs="EPSG:4326"), band_count

    return gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326"), band_count


def _laea_crs_from_bounds(minx, miny, maxx, maxy) -> PyprojCRS:
    """Local LAEA CRS (metres) for a bounding box in degrees."""
    lon_0 = (minx + maxx) / 2.0
    lat_0 = max(-89.9, min(89.9, (miny + maxy) / 2.0))
    proj4 = f"+proj=laea +lat_0={lat_0} +lon_0={lon_0} +datum=WGS84 +units=m +no_defs"
    # pyproj 2.x / 3.x / OSGeo4W builds expose different factory methods
    for factory in (
        lambda: PyprojCRS.from_user_input(proj4),
        lambda: PyprojCRS(proj4),
        lambda: PyprojCRS.from_proj4(proj4),
        lambda: PyprojCRS.from_string(proj4),
    ):
        try:
            return factory()
        except (AttributeError, TypeError, ValueError):
            continue
    raise RuntimeError(
        "Cannot construct LAEA CRS from PROJ string; upgrade pyproj or use EPSG:4326 only."
    )


def _metric_crs(gdf: gpd.GeoDataFrame) -> PyprojCRS | None:
    """Lambert azimuthal equal-area CRS centered on the data extent (metres)."""
    if gdf.crs is None:
        return None
    if not gdf.crs.is_geographic:
        return gdf.crs

    minx, miny, maxx, maxy = gdf.total_bounds
    return _laea_crs_from_bounds(minx, miny, maxx, maxy)


def _reproject_for_metric(source_gdf, target_gdf):
    """
    Reproject source and target to a shared metric CRS when inputs are geographic.

    Uses ``_metric_crs`` on the combined extent of both layers (one ``to_crs`` each).
    """
    import geopandas as gpd

    non_empty = [g for g in (source_gdf, target_gdf) if not g.empty]
    if not non_empty:
        return source_gdf, target_gdf

    minx = min(g.total_bounds[0] for g in non_empty)
    miny = min(g.total_bounds[1] for g in non_empty)
    maxx = max(g.total_bounds[2] for g in non_empty)
    maxy = max(g.total_bounds[3] for g in non_empty)
    extent = gpd.GeoDataFrame(
        geometry=[box(minx, miny, maxx, maxy)], crs=non_empty[0].crs
    )
    metric_crs = _metric_crs(extent)
    if metric_crs is None:
        return source_gdf, target_gdf
    if source_gdf.crs is not None and source_gdf.crs.equals(metric_crs):
        return source_gdf, target_gdf
    return source_gdf.to_crs(metric_crs), target_gdf.to_crs(metric_crs)


def footprint_gdf_from_raster(raster_path: str):
    """Minimal GeoDataFrame covering the raster extent (for bbox-scoped grid generation)."""
    import geopandas as gpd

    bbox = raster_bbox_wgs84(raster_path)
    return gpd.GeoDataFrame(
        geometry=[box(bbox[0], bbox[1], bbox[2], bbox[3])],
        crs="EPSG:4326",
    )


def nearest_neighbour_from_grid(raster_path: str, grid_gdf):
    """
    Assign each grid cell the band values of the nearest raster pixel center.

    Uses ``_metric_crs`` (local LAEA) when the grid is geographic.
    """
    import geopandas as gpd

    if grid_gdf.empty:
        raise ValueError("No grid cells were generated for the raster extent.")

    pixel_gdf, band_count = read_pixel_centroids(raster_path)
    if pixel_gdf.empty:
        raise ValueError("No valid raster pixels found.")

    metric_crs = _metric_crs(grid_gdf)
    if (
        metric_crs is not None
        and grid_gdf.crs is not None
        and grid_gdf.crs.is_geographic
    ):
        grid_metric = grid_gdf.to_crs(metric_crs)
        pixel_metric = pixel_gdf.to_crs(metric_crs)
    else:
        grid_metric = grid_gdf
        pixel_metric = pixel_gdf

    out = grid_gdf.copy()
    grid_pts = gpd.GeoDataFrame(
        index=grid_metric.index,
        geometry=grid_metric.geometry.centroid,
        crs=grid_metric.crs,
    )
    for bi in range(band_count):
        field = f"band_{bi + 1}"
        source_pts = gpd.GeoDataFrame(
            {field: pixel_metric[field].values},
            geometry=pixel_metric.geometry,
            crs=pixel_metric.crs,
        )
        nearest = gpd.sjoin_nearest(grid_pts, source_pts, how="left")
        out[field] = nearest[field].astype(float).round(3)

    if band_count > 0:
        out = out[out["band_1"].notna()].copy()

    return out


def accumulate_raster_pixels(
    raster_path: str,
    cell_id_fn,
    stats: str,
    desc: str = "Binning raster blocks",
):
    """
    Stream raster blocks and aggregate pixel values per cell id from ``cell_id_fn(lat, lon)``.

    Returns
    -------
    tuple
        ``(accumulator_dict, band_count)``
    """
    import rasterio
    from pyproj import Transformer
    from tqdm import tqdm

    acc = {}
    band_count = 0

    with rasterio.open(raster_path) as src:
        if src.crs is None:
            raise ValueError(
                "Raster CRS is undefined. Raster conversion requires a valid CRS."
            )
        transform = src.transform
        band_count = src.count
        transformer = None
        if src.crs.to_epsg() != 4326:
            transformer = Transformer.from_crs(src.crs, "EPSG:4326", always_xy=True)

        block_height, block_width = src.block_shapes[0]
        total_blocks = ((src.height + block_height - 1) // block_height) * (
            (src.width + block_width - 1) // block_width
        )
        for _, window in tqdm(
            src.block_windows(1),
            total=total_blocks,
            desc=desc,
            unit=" block",
        ):
            block = src.read(window=window, masked=True)
            if block.size == 0:
                continue

            block_mask = np.ma.getmaskarray(block)
            valid_pixels = ~np.all(block_mask, axis=0)
            if not np.any(valid_pixels):
                continue

            rows_local, cols_local = np.where(valid_pixels)
            rows_global = rows_local + window.row_off
            cols_global = cols_local + window.col_off

            xs, ys = rasterio.transform.xy(
                transform, rows_global, cols_global, offset="center"
            )
            if transformer is not None:
                lons, lats = transformer.transform(xs, ys)
            else:
                lons, lats = xs, ys

            for idx, (lon, lat) in enumerate(zip(lons, lats)):
                pixel_vals = block[:, rows_local[idx], cols_local[idx]]
                cell_id = cell_id_fn(lat, lon)
                if cell_id is None:
                    continue
                if cell_id not in acc:
                    acc[cell_id] = init_dggs_band_accumulator(band_count, stats)

                clean_vals = []
                for value in pixel_vals:
                    if np.ma.is_masked(value):
                        clean_vals.append(None)
                    else:
                        clean_vals.append(value)
                update_dggs_band_accumulator(acc[cell_id], clean_vals, stats)

    return acc, band_count


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
        distance = geod.inv(lon1, lat1, lon2, lat2)[
            2
        ]  # [2] gives the distance in meters
        if distance < shortest_distance:
            shortest_distance = distance

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


def dggal_to_geo(dggs_type: str, zone_id: str, options: dict = {}):
    """
    Convert a DGGAL ZoneID to a Shapely Polygon.

    Args:
        dggs_type (str): The DGGS type e.g 'gnosis','isea4r','isea9r','isea3h','isea7h','isea7h_z7',
            'ivea4r','ivea9r','ivea3h','ivea7h','ivea7h_z7','rtea4r','rtea9r','rtea3h','rtea7h','rtea7h_z7','healpix','rhealpix'
        zone_id (str): The zone identifier
        options (dict): Options for geometry generation
            - 'crs' (str | None): "5x6" or "isea" to select alternate CRS; default WGS84

    Returns:
        shapely.geometry.Polygon: Shapely Polygon or None if invalid zone
    """
    # Validate DGGS type
    dggs_type = validate_dggal_type(dggs_type)
    # Create the appropriate DGGS instance
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()

    # Get zone
    zone = dggrs.getZoneFromTextID(zone_id)
    if zone == nullZone:
        return None

    # Get CRS option
    crsOption = options.get("crs") if options else None
    crs = 0

    if crsOption:
        if crsOption == "5x6":
            crs = CRS(ogc, 153456)
        elif crsOption == "isea":
            crs = CRS(ogc, 1534)

    # Get vertices based on CRS
    if not crs or crs == CRS(ogc, 84) or crs == CRS(epsg, 4326):
        vertices = dggrs.getZoneRefinedWGS84Vertices(zone, 0)
        if vertices and vertices.count:
            count = vertices.count
            # Create list of (lon, lat) coordinates
            coords = [
                (float(vertices[i].lon), float(vertices[i].lat)) for i in range(count)
            ]
            # Close the polygon by adding the first point at the end
            coords.append(coords[0])
            return shapely.geometry.Polygon(coords)
    else:
        vertices = dggrs.getZoneRefinedCRSVertices(zone, crs, 0)
        if vertices and vertices.count:
            count = vertices.count
            # Create list of (x, y) coordinates
            coords = [
                (float(vertices[i].x), float(vertices[i].y)) for i in range(count)
            ]
            # Close the polygon by adding the first point at the end
            coords.append(coords[0])
            return shapely.geometry.Polygon(coords)

    # Return None if no valid vertices found
    return None


def dggal_dggsjsonfile2geojson(input_file, output_file=None, options: dict = {}):
    exitCode = 1
    try:
        # Check if input is a URL
        parsed_url = urlparse(input_file)
        is_url_input = all([parsed_url.scheme, parsed_url.netloc])

        # Read the JSON file from URL or local file
        if is_url_input:
            try:
                response = requests.get(input_file)
                response.raise_for_status()
                dggal_json = response.json()
            except requests.RequestException as e:
                print(f"Failure to download file from URL {input_file}: {str(e)}")
                return exitCode
        else:
            # Read from local file
            with open(input_file, "r", encoding="utf-8") as f:
                dggal_json = json.load(f)

        if dggal_json:
            # Extract options
            centroids = options.get("centroids") if options else False
            crsOption = options.get("crs") if options else None
            crs = None

            # Convert CRS option string to CRS object
            if crsOption:
                if crsOption == "5x6":
                    crs = CRS(ogc, 153456)
                elif crsOption == "isea":
                    crs = CRS(ogc, 1534)

            # Use dggal_dggsjson2geojson to convert
            geojson_result = dggal_dggsjson2geojson(
                dggal_json, crs=crs, centroids=centroids
            )

            if geojson_result:
                # Generate output filename: <input file name>.geojson
                if is_url_input:
                    # Extract filename from URL path
                    url_path = parsed_url.path.rstrip("/")
                    if url_path:
                        input_basename = os.path.splitext(os.path.basename(url_path))[0]
                        # If no filename in URL path, use domain name
                        if not input_basename:
                            input_basename = (
                                parsed_url.netloc.replace(".", "_")
                                if parsed_url.netloc
                                else "output"
                            )
                    else:
                        # No path in URL, use domain name
                        input_basename = (
                            parsed_url.netloc.replace(".", "_")
                            if parsed_url.netloc
                            else "output"
                        )
                else:
                    # Local file: use the basename without extension
                    input_basename = os.path.splitext(os.path.basename(input_file))[0]

                if output_file is None:
                    output_file = os.path.join(os.getcwd(), f"{input_basename}.geojson")

                # Write GeoJSON to file
                with open(output_file, "w", encoding="utf-8") as f:
                    json.dump(geojson_result, f, indent=2)

                exitCode = 0
            else:
                print(f"Failed to convert DGGS-JSON to GeoJSON")
        else:
            print(f"Failure to parse DGGS-JSON file {input_file}")
    except FileNotFoundError:
        print(f"Failure to open file {input_file}")
    except json.JSONDecodeError as e:
        print(f"Failure to parse DGGS-JSON file {input_file}: {str(e)}")
    except Exception as e:
        print(f"Error processing file {input_file}: {str(e)}")

    return exitCode


def dggal_generatezonegeometry(dggrs, zone, crs, id, centroids: bool, fc: bool):
    coordinates = []
    if not crs or crs == CRS(ogc, 84) or crs == CRS(epsg, 4326):
        if centroids:
            centroid = dggrs.getZoneWGS84Centroid(zone, centroid)
            coordinates.append(centroid.lon.value)
            coordinates.append(centroid.lat.value)
        else:
            vertices = dggrs.getZoneRefinedWGS84Vertices(zone, 0)
            if vertices:
                contour = []
                for v in vertices:
                    contour.append([v.lon.value, v.lat.value])
                contour.append([vertices[0].lon.value, vertices[0].lat.value])
                coordinates.append(contour)
    else:
        if centroids:
            centroid = dggrs.getZoneCRSCentroid(zone, crs, centroid)
            coordinates.append(centroid.x)
            coordinates.append(centroid.y)
        else:
            vertices = dggrs.getZoneRefinedCRSVertices(zone, crs, 0)
            if vertices:
                count = vertices.count
                contour = []
                for v in vertices:
                    contour.append([v.x, v.y])
                contour.append([vertices[0].x, vertices[0].y])
                coordinates.append(contour)
    geometry = {"type": "Point" if centroids else "Polygon", "coordinates": coordinates}
    return geometry


def dggal_generatezonefeature(dggrs, zone, crs, id, centroids: bool, fc: bool, props):
    zoneID = dggrs.getZoneTextID(zone)

    properties = {"zoneID": f"{zoneID}"}
    if props:
        for key, v in props.items():
            properties[key] = v

    features = {
        "type": "Feature",
        "id": id if id is not None else zoneID,
        "geometry": dggal_generatezonegeometry(dggrs, zone, crs, id, centroids, fc),
        "properties": properties,
    }
    return features


def dggal_dggsjson2geojson(dggal_json, crs: PyprojCRS = None, centroids: bool = False):
    result = None
    if dggal_json is not None:
        dggrsClass = None
        dggrsID = getLastDirectory(dggal_json["dggrs"])

        # We could use globals()['GNOSISGlobalGrid'] to be more generic, but here we limit to DGGRSs we know
        if not strnicmp(dggrsID, "GNOSIS", 6):
            dggrsClass = GNOSISGlobalGrid
        elif not strnicmp(dggrsID, "ISEA4R", 6):
            dggrsClass = ISEA4R
        elif not strnicmp(dggrsID, "ISEA9R", 6):
            dggrsClass = ISEA9R
        elif not strnicmp(dggrsID, "ISEA3H", 6):
            dggrsClass = ISEA3H
        elif not strnicmp(dggrsID, "ISEA7H", 6):
            dggrsClass = ISEA7H
        elif not strnicmp(dggrsID, "IVEA4R", 6):
            dggrsClass = IVEA4R
        elif not strnicmp(dggrsID, "IVEA9R", 6):
            dggrsClass = IVEA9R
        elif not strnicmp(dggrsID, "IVEA3H", 6):
            dggrsClass = IVEA3H
        elif not strnicmp(dggrsID, "IVEA7H", 6):
            dggrsClass = IVEA7H
        elif not strnicmp(dggrsID, "RTEA4R", 6):
            dggrsClass = RTEA4R
        elif not strnicmp(dggrsID, "RTEA9R", 6):
            dggrsClass = RTEA9R
        elif not strnicmp(dggrsID, "RTEA3H", 6):
            dggrsClass = RTEA3H
        elif not strnicmp(dggrsID, "RTEA7H", 6):
            dggrsClass = RTEA7H
        elif not strnicmp(dggrsID, "HEALPix", 7):
            dggrsClass = HEALPix
        elif not strnicmp(dggrsID, "rHEALPix", 8):
            dggrsClass = rHEALPix

        if dggrsClass:
            zoneID = dggal_json["zoneId"]
            dggrs = dggrsClass()
            zone = dggrs.getZoneFromTextID(zoneID)

            if zone != nullZone:
                depths = dggal_json["depths"]
                if depths:
                    maxDepth = -1

                    for d in range(len(depths)):
                        depth = depths[d]
                        if depth > maxDepth:
                            maxDepth = depth
                            break
                    if d < len(depths):
                        depth = maxDepth
                        subZones = dggrs.getSubZones(zone, depth)
                        if subZones:
                            i = 0
                            values = dggal_json["values"]
                            features = []
                            for z in subZones:
                                props = {}

                                # NOTE: We should eventually try to support __iter__ on containers
                                #       for key, depths in dggal_json.values.items():
                                for key, vDepths in values.items():
                                    if key and vDepths and len(vDepths) > d:
                                        djDepth = vDepths[d]
                                        data = djDepth["data"]
                                        props[key] = data[i]

                                features.append(
                                    dggal_generatezonefeature(
                                        dggrs, z, crs, i + 1, centroids, True, props
                                    )
                                )
                                i += 1
                            result = {"type": "FeatureCollection", "features": features}
    return result


def characteristic_length_scale(
    cell_area: float, unit: str = "m"
):  # cell_area is in m2
    """
    Compute the Characteristic Length Scale (CLS) from the cell area.

    Parameters
    ----------
    cell_area : float
        Cell area in m²

    Returns
    -------
    float
        Characteristic length scale in the m
    """
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")
    # average cell area (always in m²)
    # cell_area = 4 * math.pi * AUTHALIC_RADIUS**2 / (N_cells)
    term = math.sqrt(cell_area / math.pi) / (2.0 * AUTHALIC_RADIUS)

    # Handle case where term > 1 (can happen with very small cell areas)
    # This indicates the cell is smaller than what can be represented on the sphere
    if term > 1:
        # Return a very small value as the characteristic length scale
        # This is the maximum possible CLS (half the circumference of the Earth)
        cls = 2.0 * math.pi * AUTHALIC_RADIUS
    else:
        cls = 4.0 * AUTHALIC_RADIUS * math.asin(term)
    if unit == "km":
        cls = cls / (10**3)
    return cls


def characteristic_length_scale_from_num_cells(N_cells: int, unit: str = "m"):
    """
    Compute the Characteristic Length Scale (CLS) from the number of cells.

    Parameters
    ----------
    N_cells : int
        Number of cells at the given resolution.

    Returns
    -------
    float
        Characteristic length scale in m

    """
    # Validate unit parameter
    unit = unit.strip().lower()
    if unit not in {"m", "km"}:
        raise ValueError("unit must be one of {'m','km'}")

    # average cell area (always in m²)
    cell_area = 4 * math.pi * AUTHALIC_RADIUS**2 / (N_cells)
    term = math.sqrt(cell_area / math.pi) / (2.0 * AUTHALIC_RADIUS)

    # Handle case where term > 1 (can happen with very large numbers of cells)
    # This indicates the cell is smaller than what can be represented on the sphere
    if term > 1:
        # Return a very small value as the characteristic length scale
        # This is the maximum possible CLS (half the circumference of the Earth)
        cls = 2.0 * math.pi * AUTHALIC_RADIUS
    else:
        cls = 4.0 * AUTHALIC_RADIUS * math.asin(term)

    if unit == "km":
        cls = cls / (10**3)
    return cls


def dggrid_num_edges(dggs_type: str) -> int:
    """Nominal edge count per DGGRID DGGS type for geodesic_dggs_metrics."""
    dggs_type = validate_dggrid_type(dggs_type)
    if dggs_type in ("ISEA4T", "FULLER4T", "ISEA4D", "FULLER4D"):
        return 4
    return 6


def get_area_perimeter_from_lambert(geom):
    """Area from cell's lambert azimutal projection"""
    perimeter = np.nan
    area = np.nan
    try:
        if (-180 <= geom.centroid.x <= 180) and (-90 <= geom.centroid.y <= 90):
            proj_str = f"+proj=laea +lat_0={geom.centroid.y} +lon_0={geom.centroid.x}"
            project = pyproj.Transformer.from_crs(
                "EPSG:4236", proj_str, always_xy=True
            ).transform
            perimeter = transform(project, geom).length
            area = transform(project, geom).area
        # else:
        #     print(f'invalid centroid {geom.centroid}')
    except Exception as ex:
        # print(f'invalid centroid Exception')
        print(ex)
    # return (area, perimeter)
    return pd.Series([area, perimeter])


def convexhull_from_lambert(geom):
    """
    Compute the convex hull of a geometry using Lambert Azimuthal Equal Area projection.

    The geometry is projected to LAEA centered at its centroid, the convex hull is computed
    in the projected space, and then transformed back to WGS84 (EPSG:4326).

    Args:
        geom: Shapely geometry (Polygon, MultiPolygon, Point, etc.)

    Returns:
        Shapely Polygon: Convex hull in WGS84 coordinates, or None if invalid
    """
    try:
        if geom is None or geom.is_empty:
            return None

        # Validate centroid coordinates
        if not (-180 <= geom.centroid.x <= 180) or not (-90 <= geom.centroid.y <= 90):
            # print(f'invalid centroid {geom.centroid}')
            return None

        # Create LAEA projection centered at geometry centroid
        proj_str = f"+proj=laea +lat_0={geom.centroid.y} +lon_0={geom.centroid.x}"

        # Create transformers: WGS84 -> LAEA and LAEA -> WGS84
        to_lambert = pyproj.Transformer.from_crs(
            "EPSG:4326", proj_str, always_xy=True
        ).transform
        from_lambert = pyproj.Transformer.from_crs(
            proj_str, "EPSG:4326", always_xy=True
        ).transform

        # Transform geometry to LAEA
        geom_lambert = transform(to_lambert, geom)

        # Compute convex hull in LAEA space
        hull_lambert = geom_lambert.convex_hull

        # Transform convex hull back to WGS84
        hull_wgs84 = transform(from_lambert, hull_lambert)

        return hull_wgs84

    except Exception as ex:
        print(f"Error computing convex hull from Lambert: {ex}")
        return None


def compactness_calculation_np(df, spherical=False):
    np_area = df["area"].values
    np_peri = df["perimeter"].values

    if not spherical:
        top = 4 * np.pi * np_area
        bottom = np.square(np_peri)
        c_orig = top / bottom

        return c_orig
    else:
        t1 = 4 * np.pi * np_area
        t2 = np.square(np_area)
        t3 = np.square(6378137.0)
        top = np.sqrt(t1 - t2 / t3)
        bottom = np_peri
        zsc_np = top / bottom

        # zsc_orig = df.apply(lambda row: math.sqrt(  4*math.pi*row['area'] - math.pow(row['area'],2) / math.pow(6378137,2)  )  /   row['perimeter'], axis=1).values
        # diffs = np.sum(zsc_np - zsc_orig)
        # ma = max_error(zsc_orig, zsc_np)
        # mae = mean_absolute_error(zsc_orig, zsc_np)
        # print(f"ma={ma} mae={mae} diffs={diffs}")

        return zsc_np


def zsc_calculation_np(df, spherical=False):
    np_area = df["area"].values
    np_peri = df["perimeter"].values

    if not spherical:
        top = 4 * np.pi * np_area
        bottom = np.square(np_peri)
        c_orig = top / bottom

        return c_orig
    else:
        t1 = 4 * np.pi * np_area
        t2 = np.square(np_area)
        t3 = np.square(6378137.0)
        top = np.sqrt(t1 - t2 / t3)
        bottom = np_peri
        zsc_np = top / bottom

        # zsc_orig = df.apply(lambda row: math.sqrt(  4*math.pi*row['area'] - math.pow(row['area'],2) / math.pow(6378137,2)  )  /   row['perimeter'], axis=1).values
        # diffs = np.sum(zsc_np - zsc_orig)
        # ma = max_error(zsc_orig, zsc_np)
        # mae = mean_absolute_error(zsc_orig, zsc_np)
        # print(f"ma={ma} mae={mae} diffs={diffs}")

        return zsc_np


def get_cells_area(gdf, crs):
    """Get cells area for crs"""

    if crs == "LAEA":
        # gdf['area'],gdf['perimeter'] = zip(*gdf['geometry'].apply(get_area_perimeter_from_lambert))
        gdf[["area", "perimeter"]] = gdf["geometry"].apply(
            get_area_perimeter_from_lambert
        )

    else:
        gdf = gdf.to_crs(crs)
        gdf["area"] = gdf["geometry"].area
    return gdf


def get_a5_resolution_from_area(meters_squared, relative_depth=8):
    min_res = DGGS_TYPES["a5"]["min_res"]
    max_res = (
        DGGS_TYPES["a5"]["max_res"] + relative_depth
    )  # https://dggal.org/docs/html/dggal/Classes/DGGRS/Methods/get64KDepth.html
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    cell_count = 12
    for res in range(min_res, max_res + 1):
        # cell_count = a5.get_num_cells(res)
        if res >= 1:
            cell_count = 15 * 4 ** (res)
        if cell_count >= target_cell_count:
            return res
    return res


def get_a5_resolution_from_scale_denominator(
    scale_denominator, relative_depth=8, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_a5_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_h3_resolution_from_area(meters_squared, relative_depth=6):
    min_res = DGGS_TYPES["h3"]["min_res"]
    max_res = DGGS_TYPES["h3"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    for res in range(min_res, max_res + 1):
        cell_count = 2 + 120 * 7 ** (res)
        if cell_count >= target_cell_count:
            return res
    return res


def get_h3_resolution_from_scale_denominator(
    scale_denominator, relative_depth=6, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_h3_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_s2_resolution_from_area(meters_squared, relative_depth=8):
    min_res = DGGS_TYPES["s2"]["min_res"]
    max_res = DGGS_TYPES["s2"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    for res in range(min_res, max_res + 1):
        cell_count = 4 ** (res)
        if cell_count >= target_cell_count:
            return res
    return res


def get_s2_resolution_from_scale_denominator(
    scale_denominator, relative_depth=8, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_s2_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_rhealpix_resolution_from_area(meters_squared, relative_depth=5):
    min_res = DGGS_TYPES["rhealpix"]["min_res"]
    max_res = DGGS_TYPES["rhealpix"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    for res in range(min_res, max_res + 1):
        cell_count = 6 * 9 ** (res)
        if cell_count >= target_cell_count:
            return res
    return res


def get_rhealpix_resolution_from_scale_denominator(
    scale_denominator, relative_depth=5, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_rhealpix_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_isea4t_resolution_from_area(meters_squared, relative_depth=8):
    min_res = DGGS_TYPES["isea4t"]["min_res"]
    max_res = DGGS_TYPES["isea4t"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    for res in range(min_res, max_res + 1):
        cell_count = 20 * (4**res)
        if cell_count >= target_cell_count:
            return res
    return res


def get_isea4t_resolution_from_scale_denominator(
    scale_denominator, relative_depth=8, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_isea4t_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_isea3h_resolution_from_area(meters_squared, relative_depth=10):
    min_res = DGGS_TYPES["isea3h"]["min_res"]
    max_res = DGGS_TYPES["isea3h"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    for res in range(min_res, max_res + 1):
        cell_count = 10 * (3**res) + 2
        if cell_count >= target_cell_count:
            return res
    return res


def get_isea3h_resolution_from_scale_denominator(
    scale_denominator, relative_depth=10, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_isea3h_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_geohash_resolution_from_area(meters_squared, relative_depth=3):
    min_res = DGGS_TYPES["geohash"]["min_res"]
    max_res = DGGS_TYPES["geohash"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0
    for res in range(min_res, max_res + 1):
        cell_count = 32**res
        if cell_count >= target_cell_count:
            return res
    return res


def get_qtm_resolution_from_area(meters_squared, relative_depth=8):
    min_res = DGGS_TYPES["qtm"]["min_res"]
    max_res = DGGS_TYPES["qtm"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 1
    for res in range(min_res, max_res + 1):
        cell_count = 8 * 4 ** (res - 1)
        if cell_count >= target_cell_count:
            return res
    return res


def get_qtm_resolution_from_scale_denominator(
    scale_denominator, relative_depth=8, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_qtm_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_geohash_resolution_from_scale_denominator(
    scale_denominator, relative_depth=3, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        1,
        get_geohash_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_georef_resolution_from_area(meters_squared, relative_depth=2):
    min_res = DGGS_TYPES["georef"]["min_res"]
    max_res = DGGS_TYPES["georef"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0

    for res in range(min_res, max_res + 1):
        grid_size_deg = GEOREF_RESOLUTION_DEGREES.get(
            res
        )  # Number of cells across longitude and latitude
        num_lon = int(round(360.0 / grid_size_deg))
        num_lat = int(round(180.0 / grid_size_deg))
        cell_count = num_lon * num_lat
        if cell_count >= target_cell_count:
            return res
    return res


def get_georef_resolution_from_scale_denominator(
    scale_denominator, relative_depth=2, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_georef_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_tilecode_resolution_from_area(meters_squared, relative_depth=8):
    min_res = DGGS_TYPES["tilecode"]["min_res"]
    max_res = DGGS_TYPES["tilecode"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0

    for res in range(min_res, max_res + 1):
        cell_count = 4**res
        if cell_count >= target_cell_count:
            return res
    return res


def get_tilecode_resolution_from_scale_denominator(
    scale_denominator, relative_depth=8, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        0,
        get_tilecode_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )


def get_maidenhead_resolution_from_area(meters_squared, relative_depth=2):
    min_res = DGGS_TYPES["maidenhead"]["min_res"]
    max_res = DGGS_TYPES["maidenhead"]["max_res"] + relative_depth
    target_cell_count = AUTHALIC_AREA / meters_squared
    res = 0

    for res in range(min_res, max_res + 1):
        cell_count = maidenhead.num_cells(res)
        if cell_count >= target_cell_count:
            return res
    return res


def get_maidenhead_resolution_from_scale_denominator(
    scale_denominator, relative_depth=2, mm_per_pixel=0.28
):
    display_meters_per_pixel = (
        mm_per_pixel / 1000.0 if mm_per_pixel else STANDARD_METERS_PER_PIXEL
    )
    physical_meters_per_cell = scale_denominator * display_meters_per_pixel
    return max(
        1,
        get_maidenhead_resolution_from_area(
            physical_meters_per_cell * physical_meters_per_cell, relative_depth
        )
        - relative_depth,
    )

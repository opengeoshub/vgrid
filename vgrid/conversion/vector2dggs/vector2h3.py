import sys
import argparse
from tqdm import tqdm
from shapely.geometry import Polygon, box
import h3
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from pyproj import Geod
from vgrid.utils.geometry import calculate_point_distances
import pandas as pd
import geopandas as gpd
from vgrid.utils.geometry import (
    densify_line,
    geodesic_buffer,
    geodesic_distance,
    check_predicate,
)

geod = Geod(ellps="WGS84")

def validate_h3_resolution(resolution):
    """
    Validate that H3 resolution is in the valid range [0..15].

    Args:
        resolution: Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [0..15]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 0 or resolution > 15:
        raise ValueError(f"Resolution must be in range [0..15], got {resolution}")

    return resolution


def get_nearest_h3_resolution(points):
    """
    Find the first H3 resolution where avg_edge_length < shortest_distance.
    Uses h3.average_edge_length(resolution) to compare with the calculated shortest distance.

    Args:
        points: Shapely MultiPoint geometry

    Returns:
        int: First H3 resolution where avg_edge_length < shortest_distance (0-15)
    """
    # Calculate the shortest distance between points
    shortest_distance = calculate_point_distances(points)

    # If only one point or no distance, return a default resolution
    if shortest_distance == 0:
        return 0  # Default resolution for single points

    # Check resolutions from 0 to 15 (lowest to highest)
    for resolution in range(16):
        # Get average edge length for this H3 resolution in meters
        avg_edge_length = h3.average_hexagon_edge_length(res=resolution, unit="m")

        # Return the first resolution where avg_edge_length < shortest_distance
        if avg_edge_length <= shortest_distance:
            return resolution

    # If no resolution found where avg_edge_length < shortest_distance, return the highest resolution
    return 15


# Function to generate grid for Point
def point2h3(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a point to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        point: Shapely Point geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    h3_features = []
    # Convert point to the seed cell
    h3_id = h3.latlng_to_cell(point.y, point.x, resolution)

    cell_boundary = h3.cell_to_boundary(h3_id)
    # Wrap and filter the boundary
    filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
    # Reverse lat/lon to lon/lat for GeoJSON compatibility
    reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
    cell_polygon = Polygon(reversed_boundary)
    if cell_polygon:
        num_edges = 6
        if h3.is_pentagon(h3_id):
            num_edges = 5
        h3_feature = geodesic_dggs_to_feature(
            "h3", h3_id, resolution, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            h3_feature["properties"].update(feature_properties)
        h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def polyline2h3_new(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polyline to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        feature: Shapely LineString or MultiLineString geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    h3_features = []
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        avg_edge_length = h3.average_hexagon_edge_length(resolution, unit="m")

        # Get the centroid of the feature to use for latitude-dependent calculation
        feature_centroid = polyline.centroid
        feature_lat = feature_centroid.y
        feature_lon = feature_centroid.x

        # Calculate segment length in degrees using pyproj Geod
        segment_length_degrees, _ = geodesic_distance(
            feature_lat, feature_lon, avg_edge_length
        )

        densified_geometry = densify_line(polyline, segment_length_degrees)

        # Extract coordinates from the densified LineString
        vertices = list(densified_geometry.coords)  # (lon, lat) format

        if len(vertices) < 2:
            continue

        # Convert to (lat, lon) format for H3
        h3_vertices = [(lat, lon) for lon, lat in vertices]  # (lat, lon)

        # h3_cells = set()
        h3_cells = []
        for lat, lon in h3_vertices:
            h3_cell = h3.latlng_to_cell(lat, lon, resolution)
            h3_cells.append(h3_cell)

        for h3_cell in h3_cells:
            cell_boundary = h3.cell_to_boundary(h3_cell)
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            cell_resolution = h3.get_resolution(h3_cell)
            num_edges = 6
            if h3.is_pentagon(h3_cell):
                num_edges = 5
            h3_feature = geodesic_dggs_to_feature(
                "h3", h3_cell, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                h3_feature["properties"].update(feature_properties)
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def polyline2h3(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polyline to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        feature: Shapely LineString or MultiLineString geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    h3_features = []
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        bbox = box(*polyline.bounds)
        distance = h3.average_hexagon_edge_length(resolution, unit="m") * 2
        bbox_buffer = geodesic_buffer(bbox, distance)
        bbox_buffer_cells = h3.geo_to_cells(bbox_buffer, resolution)
        if compact:
            bbox_buffer_cells = h3.compact_cells(bbox_buffer_cells)

        for bbox_buffer_cell in bbox_buffer_cells:
            cell_boundary = h3.cell_to_boundary(bbox_buffer_cell)
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)

            # Use the check_predicate function to determine if we should keep this cell
            if not check_predicate(cell_polygon, polyline, "intersects"):
                continue  # Skip non-matching cells

            cell_resolution = h3.get_resolution(bbox_buffer_cell)
            num_edges = 6
            if h3.is_pentagon(bbox_buffer_cell):
                num_edges = 5
            h3_feature = geodesic_dggs_to_feature(
                "h3", bbox_buffer_cell, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                h3_feature["properties"].update(feature_properties)
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def polygon2h3(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polygon to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        feature: Shapely Polygon or MultiPolygon geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    h3_features = []
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        bbox = box(*polygon.bounds)
        distance = h3.average_hexagon_edge_length(resolution, unit="m") * 2
        bbox_buffer = geodesic_buffer(bbox, distance)
        bbox_buffer_cells = h3.geo_to_cells(bbox_buffer, resolution)
        if compact:
            bbox_buffer_cells = h3.compact_cells(bbox_buffer_cells)

        for bbox_buffer_cell in bbox_buffer_cells:
            cell_boundary = h3.cell_to_boundary(bbox_buffer_cell)
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)

            # Use the check_predicate function to determine if we should keep this cell
            if not check_predicate(cell_polygon, polygon, predicate):
                continue  # Skip non-matching cells

            cell_resolution = h3.get_resolution(bbox_buffer_cell)
            num_edges = 6
            if h3.is_pentagon(bbox_buffer_cell):
                num_edges = 5
            h3_feature = geodesic_dggs_to_feature(
                "h3", bbox_buffer_cell, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                h3_feature["properties"].update(feature_properties)
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def geometry2h3(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert Shapely geometry objects directly to H3 grid cells without converting to GeoJSON first.

    Args:
        geometries: Single Shapely geometry or list of Shapely geometries
        resolution (int): H3 resolution [0..15]
        properties_list: List of property dictionaries (optional)
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode - for polygon only
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    resolution = validate_h3_resolution(resolution)

    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    h3_features = []

    for i, geom in enumerate(tqdm(geometries, desc="Processing features")):
        if geom is None:
            continue

        # Get properties for this geometry
        props = properties_list[i] if i < len(properties_list) else {}

        # Process based on geometry type
        if geom.geom_type == "Point":
            point_features = point2h3(
                resolution,
                geom,
                props,
                predicate,
                compact,
                topology,
                include_properties,
            )
            h3_features.extend(point_features["features"])

        elif geom.geom_type == "MultiPoint":
            for point in geom.geoms:
                point_features = point2h3(
                    resolution,
                    point,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
                h3_features.extend(point_features["features"])

        elif geom.geom_type in ["LineString", "MultiLineString"]:
            polyline_features = polyline2h3(
                resolution,
                geom,
                props,
                predicate,
                compact,
                topology,
                include_properties,
            )
            h3_features.extend(polyline_features["features"])

        elif geom.geom_type in ["Polygon", "MultiPolygon"]:
            poly_features = polygon2h3(
                resolution,
                geom,
                props,
                predicate,
                compact,
                topology,
                include_properties,
            )
            h3_features.extend(poly_features["features"])

        else:
            raise ValueError(f"Unsupported geometry type: {geom.geom_type}")

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


def dataframe2h3(
    df,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert pandas DataFrame with geometry column to H3 grid cells by converting to Shapely geometries first.

    Args:
        df (pd.DataFrame): Input DataFrame with geometry column
        resolution (int): H3 resolution [0..15]
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode - for polygon only
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    # Find geometry column
    geometry_col = None
    for col in df.columns:
        if hasattr(df[col].iloc[0], "geom_type") or hasattr(
            df[col].iloc[0], "__geo_interface__"
        ):
            geometry_col = col
            break

    if geometry_col is None:
        raise ValueError(
            "DataFrame must contain a geometry column with Shapely geometry objects"
        )

    # Extract geometries and properties from DataFrame
    geometries = []
    properties_list = []

    for idx, row in df.iterrows():
        # Get the geometry
        geom = row[geometry_col]
        if geom is not None:
            geometries.append(geom)

            # Get properties (exclude geometry column)
            properties = row.to_dict()
            if geometry_col in properties:
                del properties[geometry_col]
            properties_list.append(properties)

    # Use geometry2h3 to process the geometries
    return geometry2h3(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2h3(
    gdf,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert GeoDataFrame to H3 grid cells by converting to Shapely geometries first.

    Args:
        gdf (GeoDataFrame): Input GeoDataFrame
        resolution (int): H3 resolution [0..15]
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode - for polygon only
        topology (bool): Enable H3 topology preserving mode
        include_properties (bool): If False, do not include original feature properties

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    # Extract geometries and properties from GeoDataFrame
    geometries = []
    properties_list = []

    for idx, row in gdf.iterrows():
        # Get the geometry
        geom = row.geometry
        if geom is not None:
            geometries.append(geom)

            # Get properties (exclude geometry column)
            properties = row.to_dict()
            if "geometry" in properties:
                del properties["geometry"]
            properties_list.append(properties)

    # Use geometry2h3 to process the geometries
    return geometry2h3(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def vector2h3(
    data,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    output_format="geojson",
    output_path=None,
    include_properties=True,
    **kwargs,
):
    """
    Convert vector data to H3 grid cells from various input formats.

    This function can handle:
    - File paths (GeoJSON, Shapefile, GeoPackage, KML, GML, etc.) - converted to GeoDataFrame then Shapely geometries
    - URLs (remote files) - converted to GeoDataFrame then Shapely geometries
    - pandas DataFrames with geometry column - converted to Shapely geometries
    - GeoJSON dictionaries - converted to GeoDataFrame then Shapely geometries
    - Shapely geometry objects - processed directly
    - GeoDataFrames - converted to Shapely geometries

    Args:
        data: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
        resolution (int): H3 resolution [0..15]
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode for polygons (default: False)
        topology (bool): Enable H3 topology preserving mode (default: False)
        output_format (str): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str): Output file path (optional)
        include_properties (bool): If False, do not include original feature properties. (default: True)
        **kwargs: Additional arguments passed to geopandas read functions

    Returns:
        dict or str: Output in the specified format
    """
    # Process input data directly
    if hasattr(data, "geometry") and hasattr(data, "columns"):
        # GeoDataFrame - convert to Shapely geometries first
        result = geodataframe2h3(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        # Regular DataFrame with geometry column - convert to Shapely geometries first
        result = dataframe2h3(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        # Shapely geometry objects - process directly
        result = geometry2h3(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        # GeoJSON data - convert to GeoDataFrame first, then process
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2h3(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    
    elif isinstance(data, str):
        # File path or URL - use geopandas.read_file directly
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2h3(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")

    # Convert result to specified output format
    return convert_to_output_format(result, output_format, output_path)


def convert_to_output_format(result, output_format, output_path=None):
    """
    Convert H3 result to specified output format.

    Args:
        result (dict): GeoJSON FeatureCollection result
        output_format (str): Desired output format
        output_path (str): Output file path (optional)

    Returns:
        dict or str: Output in the specified format
    """
    # First convert GeoJSON result to GeoDataFrame
    gdf = gpd.GeoDataFrame.from_features(result["features"])

    # Set CRS to WGS84 (EPSG:4326) since H3 uses WGS84 coordinates
    gdf.set_crs(epsg=4326, inplace=True)

    if output_format.lower() == "geojson":
        if output_path:
            import json

            with open(output_path, "w") as f:
                json.dump(result, f, indent=2)
            return output_path
        else:
            return result  # Already in GeoJSON format

    elif output_format.lower() == "gpkg":
        if output_path:
            gdf.to_file(output_path, driver="GPKG")
            return output_path
        else:
            gdf.to_file("vector2h3.gpkg", driver="GPKG")
            return "vector2h3.gpkg"

    elif output_format.lower() == "parquet":
        if output_path:
            gdf.to_parquet(output_path, index=False)
            return output_path
        else:
            gdf.to_parquet("vector2h3.parquet", index=False)
            return "vector2h3.parquet"

    elif output_format.lower() == "csv":
        if output_path:
            gdf.to_csv(output_path, index=False)
            return output_path
        else:
            return gdf.to_csv(index=False)

    elif output_format.lower() == "shapefile":
        if output_path:
            gdf.to_file(output_path, driver="ESRI Shapefile")
            return output_path
        else:
            gdf.to_file("vector2h3.shp", driver="ESRI Shapefile")
            return "vector2h3.shp"

    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


def vector2h3_cli():
    """Command-line interface for vector2h3 conversion."""
    parser = argparse.ArgumentParser(description="Convert vector data to H3 grid cells")
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(16),
        metavar="[0-15]",
        help="H3 resolution [0..15] (0=coarsest, 15=finest)",
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate: intersect, within, centroid_within, largest_overlap for polygons",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable H3 compact mode for polygons",
    )
    parser.add_argument(
        "-t", "--topology", action="store_true", help="Enable topology preserving mode"
    )
    parser.add_argument(
        "-np",
        "-no-props",
        dest="include_properties",
        action="store_false",
        help="Do not include original feature properties.",
    )
    parser.add_argument(
        "-f",
        "--format",
        default="geojson",
        choices=["geojson", "gpkg", "parquet", "csv", "shapefile"],
        help="Output format (default: geojson)",
    )
    parser.add_argument("-o", "--output", help="Output file path (optional)")

    args = parser.parse_args()

    # Validate resolution if provided
    if args.resolution is not None:
        try:
            args.resolution = validate_h3_resolution(args.resolution)
        except (ValueError, TypeError) as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    # Handle input (no stdin support)
    data = args.input

    # Handle output path
    output_path = args.output

    # If output path is not specified but format requires a file, generate default name
    if not output_path and args.format in [
        "geojson",
        "gpkg",
        "parquet",
        "csv",
        "shapefile",
    ]:
        # Generate default filename based on format
        extensions = {
            "geojson": ".geojson",
            "gpkg": ".gpkg",
            "parquet": ".parquet",
            "csv": ".csv",
            "shapefile": ".shp",
        }
        output_path = f"vector2h3{extensions.get(args.format, '')}"

    try:
        vector2h3(
            data,
            args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            output_format=args.format,
            output_path=output_path,
            include_properties=args.include_properties,
        )

        if output_path:
            print(f"Output saved to {output_path}")

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2h3_cli()

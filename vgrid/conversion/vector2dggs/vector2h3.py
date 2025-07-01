"""Convert vector data to H3 grid cells."""

import sys
import argparse
from tqdm import tqdm
from pyproj import Geod
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, box,MultiPoint, MultiLineString, MultiPolygon
import h3
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.geometry import (
    geodesic_buffer,
    check_predicate,
    shortest_point_distance,
    shortest_polyline_distance,
    shortest_polygon_distance,
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


# Function to generate grid for Point
def point2h3(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    all_points=None,  # New parameter for topology preservation
):
    """
    Convert a point to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        point: Shapely Point geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode - ensures disjoint points have disjoint H3 cells
        include_properties (bool): If False, do not include original feature properties
        all_points: List of all points for topology preservation (required when topology=True)

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_points is None:
            raise ValueError("all_points parameter is required when topology=True")
        
        # Calculate the shortest distance between all points
        shortest_distance = shortest_point_distance(all_points)
        
        # Find resolution where H3 cell size is smaller than shortest distance
        # This ensures disjoint points have disjoint H3 cells
        if shortest_distance > 0:
            for res in range(16):
                avg_edge_length = h3.average_hexagon_edge_length(res=res, unit="m")
                # Use a factor to ensure sufficient separation (hexagon diameter is ~2x edge length)
                hexagon_diameter = avg_edge_length*2 
                if hexagon_diameter < shortest_distance:
                    resolution = res
                    break
            else:
                # If no resolution found, use the highest resolution
                resolution = 15
        else:
            # Single point or no distance, use provided resolution
            pass

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

def polyline2h3(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    all_polylines=None,  # New parameter for topology preservation
):
    """
    Convert a polyline to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        feature: Shapely LineString or MultiLineString geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode - ensures disjoint polylines have disjoint H3 cells
        include_properties (bool): If False, do not include original feature properties
        all_polylines: List of all polylines for topology preservation (required when topology=True)

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_polylines is None:
            raise ValueError("all_polylines parameter is required when topology=True")
        
        # Calculate the shortest distance between all polylines
        shortest_distance = shortest_polyline_distance(all_polylines)
        
        # Find resolution where H3 cell size is smaller than shortest distance
        # This ensures disjoint polylines have disjoint H3 cells
        if shortest_distance > 0:
            for res in range(16):
                avg_edge_length = h3.average_hexagon_edge_length(res=res, unit="m")
                # Use a factor to ensure sufficient separation (hexagon diameter is ~2x edge length)
                hexagon_diameter = avg_edge_length * 4 # in case there are 2 cells representing the same line segment
                if hexagon_diameter < shortest_distance:
                    resolution = res
                    break
            else:
                # If no resolution found, use the highest resolution
                resolution = 15
        else:
            # Single polyline or no distance, use provided resolution
            pass

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
    all_polygons=None,  # New parameter for topology preservation
):
    """
    Convert a polygon to H3 grid cells.

    Args:
        resolution (int): H3 resolution [0..15]
        feature: Shapely Polygon or MultiPolygon geometry
        feature_properties (dict): Properties to add to the H3 features
        predicate (str or int): Spatial predicate to apply (see check_predicate function)
        compact (bool): Enable H3 compact mode
        topology (bool): Enable H3 topology preserving mode - ensures disjoint polygons have disjoint H3 cells
        include_properties (bool): If False, do not include original feature properties
        all_polygons: List of all polygons for topology preservation (required when topology=True)

    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_polygons is None:
            raise ValueError("all_polygons parameter is required when topology=True")
        
       
        # Calculate the shortest distance between all polygons
        shortest_distance = shortest_polygon_distance(all_polygons)
        
        # Find resolution where H3 cell size is smaller than shortest distance
        # This ensures disjoint polygons have disjoint H3 cells
        if shortest_distance > 0:
            for res in range(16):
                avg_edge_length = h3.average_hexagon_edge_length(res=res, unit="m")
                # Use a factor to ensure sufficient separation (hexagon diameter is ~2x edge length)
                hexagon_diameter = avg_edge_length * 4
                if hexagon_diameter < shortest_distance:
                    resolution = res
                    break
            else:
                # If no resolution found, use the highest resolution
                resolution = 15
        else:
            # Single polygon or no distance, use provided resolution
            pass

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

    # Collect all points, polylines, and polygons for topology preservation if needed
    all_points = None
    all_polylines = None
    all_polygons = None
    if topology:
        points_list = []
        polylines_list = []
        polygons_list = []
        for i, geom in enumerate(geometries):
            if geom is None:
                continue
            if geom.geom_type == "Point":
                points_list.append(geom)
            elif geom.geom_type == "MultiPoint":
                points_list.extend(list(geom.geoms))
            elif geom.geom_type == "LineString":
                polylines_list.append(geom)
            elif geom.geom_type == "MultiLineString":
                polylines_list.extend(list(geom.geoms))
            elif geom.geom_type == "Polygon":
                polygons_list.append(geom)
            elif geom.geom_type == "MultiPolygon":
                polygons_list.extend(list(geom.geoms))
        if points_list:
            all_points = MultiPoint(points_list)
        if polylines_list:
            all_polylines = MultiLineString(polylines_list)
        if polygons_list:
            all_polygons = MultiPolygon(polygons_list)

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
                all_points,  # Pass all points for topology preservation
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
                    all_points,  # Pass all points for topology preservation
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
                all_polylines,  # Pass all polylines for topology preservation
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
                all_polygons=all_polygons,  # Pass all polygons for topology preservation
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

    for _, row in df.iterrows():
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

    for _, row in gdf.iterrows():
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
    # Check if result has features
    if not result or "features" not in result or not result["features"]:
        print("Warning: No features found in result. This may happen when:")
        print("  - Using 'within' predicate with coarse resolution (cells too large)")
        print("  - Using 'largest_overlap' predicate with no cells having >50% overlap")
        print("  - Input geometry is invalid or empty")
        print("Suggestions:")
        print("  - Try a finer resolution (higher number)")
        print("  - Use 'intersect' or 'centroid_within' predicate instead")
        print("  - Check that input geometry is valid")
        raise ValueError("No features found in result")
    
    # First convert GeoJSON result to GeoDataFrame
    try:
        gdf = gpd.GeoDataFrame.from_features(result["features"])

        # Set CRS to WGS84 (EPSG:4326) since H3 uses WGS84 coordinates
        gdf.set_crs(epsg=4326, inplace=True)
        
        # Ensure the geometry column is set as the active geometry column
        if 'geometry' in gdf.columns:
            gdf.set_geometry('geometry', inplace=True)
        else:
            # If no geometry column found, try to identify it
            geom_cols = [col for col in gdf.columns if hasattr(gdf[col].iloc[0], 'geom_type')]
            if geom_cols:
                gdf.set_geometry(geom_cols[0], inplace=True)
            else:
                raise ValueError("No geometry column found in GeoDataFrame")
        
        # Verify the GeoDataFrame has valid geometry
        if gdf.empty:
            raise ValueError("GeoDataFrame is empty")
        
        if not gdf.geometry.is_valid.all():
            print("Warning: Some geometries are invalid")
    
    except Exception as e:
        print(f"Error creating GeoDataFrame: {str(e)}")
        print(f"Result features count: {len(result['features']) if 'features' in result else 0}")
        if 'features' in result and result['features']:
            print(f"First feature: {result['features'][0]}")
        raise

    if output_format.lower() == "geojson":
        if output_path:
            import json

            with open(output_path, "w", encoding="utf-8") as f:
                json.dump(result, f)
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
        "-t", "--topology", action="store_true", 
        help="Enable topology preserving mode ensuring disjoint features have disjoint H3 cells by automatically calculating appropriate resolution."
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

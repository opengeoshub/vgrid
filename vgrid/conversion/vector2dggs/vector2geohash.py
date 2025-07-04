import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Polygon, shape, MultiPoint, MultiLineString, MultiPolygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import geohash
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.generator.geohashgrid import (
    initial_geohashes,
    geohash_to_polygon,
    expand_geohash_bbox,
)
from vgrid.conversion.dggscompact import geohashcompact
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
    shortest_polyline_distance,
    shortest_polygon_distance,
)
from math import sqrt


def validate_geohash_resolution(resolution):
    """
    Validate that Geohash resolution is in the valid range [1..10] (1=coarsest, 10=finest).

    Args:
        resolution (int): Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [1..10]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 1 or resolution > 10:
        raise ValueError(f"Resolution must be in range [1..10], got {resolution}")

    return resolution


def point2geohash(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    all_points=None,
):
    """
    Convert a point geometry to a Geohash grid cell.

    Args:
        resolution (int): Geohash resolution [1..10]
        point (shapely.geometry.Point): Point geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for points)
        compact (bool, optional): Enable Geohash compact mode (not used for points)
        topology (bool, optional): Enable topology preserving mode - ensures disjoint points have disjoint Geohash cells
        include_properties (bool, optional): Whether to include properties in output
        all_points: List of all points for topology preservation (required when topology=True)

    Returns:
        list: List of GeoJSON feature dictionaries representing Geohash cells containing the point
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_points is None:
            raise ValueError("all_points parameter is required when topology=True")
        shortest_distance = shortest_point_distance(all_points)
        # Geohash cell size (approx): each character increases resolution, cell size shrinks by ~1/8 to 1/32
        # We'll use a rough estimate: cell size halves every character, so use a lookup or formula
        # For simplicity, use a fixed table for WGS84 (meters) for geohash length 1-10
        geohash_cell_sizes = [5000000, 1250000, 156000, 39100, 4890, 1220, 153, 38.2, 4.77, 1.19]  # meters
        for res in range(1, 11):
            cell_diameter = geohash_cell_sizes[res-1] * sqrt(2) * 2
            if cell_diameter < shortest_distance:
                resolution = res
                break
        else:
            resolution = 10
    geohash_features = []
    longitude = point.x
    latitude = point.y
    geohash_id = geohash.encode(latitude, longitude, resolution)
    bbox = geohash.bbox(geohash_id)
    if bbox:
        min_lat, min_lon = bbox["s"], bbox["w"]
        max_lat, max_lon = bbox["n"], bbox["e"]
        cell_polygon = Polygon(
            [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat],
            ]
        )
        geohash_feature = graticule_dggs_to_feature(
            "geohash", geohash_id, resolution, cell_polygon
        )
        if include_properties and feature_properties:
            geohash_feature["properties"].update(feature_properties)
        geohash_features.append(geohash_feature)
    return geohash_features


def polyline2geohash(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    all_polylines=None,
):
    """
    Convert line geometries (LineString, MultiLineString) to Geohash grid cells.

    Args:
        resolution (int): Geohash resolution [1..10]
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Line geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for lines)
        compact (bool, optional): Enable Geohash compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode - ensures disjoint polylines have disjoint Geohash cells
        include_properties (bool, optional): Whether to include properties in output
        all_polylines: List of all polylines for topology preservation (required when topology=True)

    Returns:
        list: List of GeoJSON feature dictionaries representing Geohash cells intersecting the line
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_polylines is None:
            raise ValueError("all_polylines parameter is required when topology=True")
        shortest_distance = shortest_polyline_distance(all_polylines)
        geohash_cell_sizes = [5000000, 1250000, 156000, 39100, 4890, 1220, 153, 38.2, 4.77, 1.19]  # meters
        for res in range(1, 11):
            cell_diameter = geohash_cell_sizes[res-1] * sqrt(2) * 4  # in case there are 2 cells representing the same line segment
            if cell_diameter < shortest_distance:
                resolution = res
                break
        else:
            resolution = 10
    geohash_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        intersected_geohashes = {
            gh
            for gh in initial_geohashes
            if geohash_to_polygon(gh).intersects(polyline)
        }
        geohashes_bbox = set()
        for gh in intersected_geohashes:
            expand_geohash_bbox(gh, resolution, geohashes_bbox, polyline)
        for gh in geohashes_bbox:
            cell_polygon = geohash_to_polygon(gh)
            geohash_feature = graticule_dggs_to_feature(
                "geohash", gh, resolution, cell_polygon
            )
            if include_properties and feature_properties:
                geohash_feature["properties"].update(feature_properties)
            geohash_features.append(geohash_feature)
    geohash_geojson = {"type": "FeatureCollection", "features": geohash_features}
    if compact:
        return geohashcompact(geohash_geojson)["features"]
    return geohash_features


def polygon2geohash(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    all_polygons=None,
):
    """
    Convert polygon geometries (Polygon, MultiPolygon) to Geohash grid cells.

    Args:
        resolution (int): Geohash resolution [1..10]
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable Geohash compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode - ensures disjoint polygons have disjoint Geohash cells
        include_properties (bool, optional): Whether to include properties in output
        all_polygons: List of all polygons for topology preservation (required when topology=True)

    Returns:
        list: List of GeoJSON feature dictionaries representing Geohash cells based on predicate
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_polygons is None:
            raise ValueError("all_polygons parameter is required when topology=True")
        shortest_distance = shortest_polygon_distance(all_polygons)
        geohash_cell_sizes = [5000000, 1250000, 156000, 39100, 4890, 1220, 153, 38.2, 4.77, 1.19]  # meters
        for res in range(1, 11):
            cell_diameter = geohash_cell_sizes[res-1] * sqrt(2) * 4
            if cell_diameter < shortest_distance:
                resolution = res
                break
        else:
            resolution = 10
    geohash_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        intersected_geohashes = {
            gh for gh in initial_geohashes if geohash_to_polygon(gh).intersects(polygon)
        }
        geohashes_bbox = set()
        for gh in intersected_geohashes:
            expand_geohash_bbox(gh, resolution, geohashes_bbox, polygon)
        for gh in geohashes_bbox:
            cell_polygon = geohash_to_polygon(gh)
            geohash_feature = graticule_dggs_to_feature(
                "geohash", gh, resolution, cell_polygon
            )
            cell_geom = shape(geohash_feature["geometry"])
            if not check_predicate(cell_geom, polygon, predicate):
                continue
            if include_properties and feature_properties:
                geohash_feature["properties"].update(feature_properties)
            geohash_features.append(geohash_feature)
    geohash_geojson = {"type": "FeatureCollection", "features": geohash_features}
    if compact:
        return geohashcompact(geohash_geojson)["features"]
    return geohash_features


# --- Main geometry conversion ---
def geometry2geohash(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a list of geometries to Geohash grid cells.

    Args:
        geometries (shapely.geometry.BaseGeometry or list): Single geometry or list of geometries
        resolution (int): Geohash resolution [1..10]
        properties_list (list, optional): List of property dictionaries for each geometry
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable Geohash compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with Geohash grid cells
    """
    resolution = validate_geohash_resolution(resolution)
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

    geohash_features = []
    for i, geom in enumerate(tqdm(geometries, desc="Processing features")):
        if geom is None:
            continue
        props = properties_list[i] if i < len(properties_list) else {}
        if not include_properties:
            props = {}
        if geom.geom_type == "Point":
            geohash_features.extend(
                point2geohash(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                    all_points,  # Pass all points for topology preservation
                )
            )
        elif geom.geom_type == "MultiPoint":
            for pt in geom.geoms:
                geohash_features.extend(
                    point2geohash(
                        resolution,
                        pt,
                        props,
                        predicate,
                        compact,
                        topology,
                        include_properties,
                        all_points,  # Pass all points for topology preservation
                    )
                )
        elif geom.geom_type in ("LineString", "MultiLineString"):
            geohash_features.extend(
                polyline2geohash(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                    all_polylines,  # Pass all polylines for topology preservation
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            geohash_features.extend(
                polygon2geohash(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                    all_polygons=all_polygons,  # Pass all polygons for topology preservation
                )
            )
    return {"type": "FeatureCollection", "features": geohash_features}


# --- DataFrame/GeoDataFrame conversion ---
def dataframe2geohash(
    df,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a pandas DataFrame with geometry column to Geohash grid cells.

    Args:
        df (pandas.DataFrame): DataFrame with geometry column
        resolution (int): Geohash resolution [1..10]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable Geohash compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with Geohash grid cells
    """
    geometries = []
    properties_list = []
    for idx, row in df.iterrows():
        geom = row.geometry if "geometry" in row else row["geometry"]
        if geom is not None:
            geometries.append(geom)
            props = row.to_dict()
            if "geometry" in props:
                del props["geometry"]
            properties_list.append(props)
    return geometry2geohash(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2geohash(
    gdf,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a GeoDataFrame to Geohash grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int): Geohash resolution [1..10]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable Geohash compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with Geohash grid cells
    """
    geometries = []
    properties_list = []
    for idx, row in gdf.iterrows():
        geom = row.geometry
        if geom is not None:
            geometries.append(geom)
            props = row.to_dict()
            if "geometry" in props:
                del props["geometry"]
            properties_list.append(props)
    return geometry2geohash(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


# --- Main vector2geohash function ---
def vector2geohash(
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
    Convert vector data to Geohash grid cells from various input formats.

    Args:
        data: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
        resolution (int): Geohash resolution [1..10]
        compact (bool): Enable Geohash compact mode for polygons (default: False)
        output_format (str): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str): Output file path (optional)
        include_properties (bool): If False, do not include original feature properties. (default: True)
        **kwargs: Additional arguments passed to geopandas read functions
    Returns:
        dict or str: Output in the specified format
    """

    # Process input data directly
    if hasattr(data, "geometry") and hasattr(data, "columns"):
        # GeoDataFrame
        result = geodataframe2geohash(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        # DataFrame with geometry column
        result = dataframe2geohash(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        # Shapely geometry objects
        result = geometry2geohash(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        # GeoJSON dict
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2geohash(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        # File path or URL
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2geohash(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")
    return convert_to_output_format(result, output_format, output_path)


# --- Output format conversion ---
def convert_to_output_format(result, output_format, output_path=None):
    """
    Convert Geohash result to specified output format.

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

        # Set CRS to WGS84 (EPSG:4326) since Geohash uses WGS84 coordinates
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
            gdf.to_file("vector2geohash.gpkg", driver="GPKG")
            return "vector2geohash.gpkg"

    elif output_format.lower() == "parquet":
        if output_path:
            gdf.to_parquet(output_path, index=False)
            return output_path
        else:
            gdf.to_parquet("vector2geohash.parquet", index=False)
            return "vector2geohash.parquet"

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
            gdf.to_file("vector2geohash.shp", driver="ESRI Shapefile")
            return "vector2geohash.shp"

    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


# --- CLI ---
def vector2geohash_cli():
    """
    Command-line interface for vector2geohash conversion.

    Usage:
        python vector2geohash.py -i input.shp -r 5 -c -f geojson -o output.geojson

    Arguments:
        -i, --input: Input file path or URL
        -r, --resolution: Geohash resolution [1..10]
        -c, --compact: Enable Geohash compact mode
        -p, --predicate: Spatial predicate (intersect, within, centroid_within, largest_overlap)
        -t, --topology: Enable topology preserving mode
        -np, --no-props: Do not include original feature properties
        -f, --format: Output format (geojson, gpkg, parquet, csv, shapefile)
        -o, --output: Output file path
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to Geohash grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(1, 11),
        metavar="[1-10]",
        help="Geohash resolution [1..10] (1=coarsest, 10=finest)",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable Geohash compact mode for polygons",
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate: intersect, within, centroid_within, largest_overlap for polygons",
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
    if args.resolution is not None:
        try:
            args.resolution = validate_geohash_resolution(args.resolution)
        except (ValueError, TypeError) as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    data = args.input
    output_path = args.output
    if not output_path and args.format in [
        "geojson",
        "gpkg",
        "parquet",
        "csv",
        "shapefile",
    ]:
        extensions = {
            "geojson": ".geojson",
            "gpkg": ".gpkg",
            "parquet": ".parquet",
            "csv": ".csv",
            "shapefile": ".shp",
        }
        output_path = f"vector2geohash{extensions.get(args.format, '')}"
    try:
        vector2geohash(
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
    vector2geohash_cli()

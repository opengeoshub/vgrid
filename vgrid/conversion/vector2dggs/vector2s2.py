import sys
import argparse
import json
from tqdm import tqdm
import pandas as pd
import geopandas as gpd
from vgrid.utils import s2
from vgrid.generator.s2grid import s2_cell_to_polygon
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.geometry import check_predicate


def validate_s2_resolution(resolution):
    """
    Validate that S2 resolution is in the valid range [0..28].

    Args:
        resolution: Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [0..28]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 0 or resolution > 28:
        raise ValueError(f"Resolution must be in range [0..28], got {resolution}")

    return resolution


def point2s2(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a single point geometry to S2 grid cells.

    Args:
        resolution (int): S2 resolution level [0..28]
        point (shapely.geometry.Point): Point geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for points)
        compact (bool, optional): Enable S2 compact mode (not used for points)
        topology (bool, optional): Enable topology preserving mode (not used for points)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing S2 cells containing the point

    Example:
        >>> from shapely.geometry import Point
        >>> point = Point(-122.4194, 37.7749)  # San Francisco
        >>> cells = point2s2(10, point, {"name": "SF"})
        >>> len(cells)
        1
    """
    s2_features = []
    lat_lng = s2.LatLng.from_degrees(point.y, point.x)
    cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
    cell_id = cell_id_max_res.parent(resolution)
    s2_cell = s2.Cell(cell_id)
    cell_token = s2.CellId.to_token(s2_cell.id())
    if s2_cell:
        cell_polygon = s2_cell_to_polygon(cell_id)
        resolution = cell_id.level()
        num_edges = 4
        s2_feature = geodesic_dggs_to_feature(
            "s2", cell_token, resolution, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            s2_feature["properties"].update(feature_properties)
        s2_features.append(s2_feature)
    return s2_features


def polyline2hs2(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert line geometries (LineString, MultiLineString) to S2 grid cells.

    Args:
        resolution (int): S2 resolution level [0..28]
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Line geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for lines)
        compact (bool, optional): Enable S2 compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (not used for lines)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing S2 cells intersecting the line

    Example:
        >>> from shapely.geometry import LineString
        >>> line = LineString([(-122.4194, 37.7749), (-122.4000, 37.7800)])
        >>> cells = polyline2hs2(10, line, {"name": "route"})
        >>> len(cells) > 0
        True
    """
    s2_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        min_lng, min_lat, max_lng, max_lat = polyline.bounds
        level = resolution
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng),
        )
        covering = coverer.get_covering(region)
        cell_ids = covering
        if compact:
            covering = s2.CellUnion(covering)
            covering.normalize()
            cell_ids = covering.cell_ids()
        for cell_id in cell_ids:
            cell_polygon = s2_cell_to_polygon(cell_id)
            if not cell_polygon.intersects(polyline):
                continue
            cell_token = s2.CellId.to_token(cell_id)
            cell_resolution = cell_id.level()
            num_edges = 4
            s2_feature = geodesic_dggs_to_feature(
                "s2", cell_token, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                s2_feature["properties"].update(feature_properties)
            s2_features.append(s2_feature)
    return s2_features


def polygon2hs2(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert polygon geometries (Polygon, MultiPolygon) to S2 grid cells.

    Args:
        resolution (int): S2 resolution level [0..28]
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable S2 compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (not used for polygons)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing S2 cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2hs2(10, poly, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    s2_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        min_lng, min_lat, max_lng, max_lat = polygon.bounds
        level = resolution
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng),
        )
        covering = coverer.get_covering(region)
        cell_ids = covering
        if compact:
            covering = s2.CellUnion(covering)
            covering.normalize()
            cell_ids = covering.cell_ids()
        for cell_id in cell_ids:
            cell_polygon = s2_cell_to_polygon(cell_id)
            if not check_predicate(cell_polygon, polygon, predicate):
                continue
            cell_token = s2.CellId.to_token(cell_id)
            cell_resolution = cell_id.level()
            num_edges = 4
            s2_feature = geodesic_dggs_to_feature(
                "s2", cell_token, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                s2_feature["properties"].update(feature_properties)
            s2_features.append(s2_feature)
    return s2_features


# --- Main geometry conversion ---
def geometry2s2(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a list of geometries to S2 grid cells.

    Args:
        geometries (shapely.geometry.BaseGeometry or list): Single geometry or list of geometries
        resolution (int): S2 resolution level [0..28]
        properties_list (list, optional): List of property dictionaries for each geometry
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable S2 compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with S2 grid cells

    Example:
        >>> from shapely.geometry import Point, Polygon
        >>> geoms = [Point(-122.4194, 37.7749), Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])]
        >>> props = [{"name": "point"}, {"name": "polygon"}]
        >>> result = geometry2s2(geoms, 10, props, predicate="intersect")
        >>> result["type"]
        'FeatureCollection'
    """
    resolution = validate_s2_resolution(resolution)

    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    s2_features = []
    for idx, geom in enumerate(tqdm(geometries, desc="Processing features")):
        if geom is None:
            continue
        props = (
            properties_list[idx]
            if properties_list and idx < len(properties_list)
            else {}
        )
        if not include_properties:
            props = {}
        if geom.geom_type == "Point":
            s2_features.extend(
                point2s2(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
        elif geom.geom_type == "MultiPoint":
            for pt in geom.geoms:
                s2_features.extend(
                    point2s2(
                        resolution,
                        pt,
                        props,
                        predicate,
                        compact,
                        topology,
                        include_properties,
                    )
                )
        elif geom.geom_type in ("LineString", "MultiLineString"):
            s2_features.extend(
                polyline2hs2(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            s2_features.extend(
                polygon2hs2(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
    return {"type": "FeatureCollection", "features": s2_features}


# --- DataFrame/GeoDataFrame conversion ---
def dataframe2s2(
    df,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a pandas DataFrame with geometry column to S2 grid cells.

    Args:
        df (pandas.DataFrame): DataFrame with geometry column
        resolution (int): S2 resolution level [0..28]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable S2 compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with S2 grid cells

    Example:
        >>> import pandas as pd
        >>> from shapely.geometry import Point
        >>> df = pd.DataFrame({
        ...     'geometry': [Point(-122.4194, 37.7749)],
        ...     'name': ['San Francisco']
        ... })
        >>> result = dataframe2s2(df, 10)
        >>> result["type"]
        'FeatureCollection'
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
    return geometry2s2(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2s2(
    gdf,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a GeoDataFrame to S2 grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int): S2 resolution level [0..28]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable S2 compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with S2 grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2s2(gdf, 10)
        >>> result["type"]
        'FeatureCollection'
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
    return geometry2s2(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


# --- Main vector2s2 function ---
def vector2s2(
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
    Convert vector data to S2 grid cells from various input formats.

    This is the main function that handles conversion of vector data to S2 grid cells.
    It supports multiple input formats including file paths, URLs, DataFrames, GeoDataFrames,
    GeoJSON dictionaries, and Shapely geometries.

    Args:
        data: Input data in one of the following formats:
            - File path (str): Path to vector file (shapefile, GeoJSON, etc.)
            - URL (str): URL to vector data
            - pandas.DataFrame: DataFrame with geometry column
            - geopandas.GeoDataFrame: GeoDataFrame
            - dict: GeoJSON dictionary
            - shapely.geometry.BaseGeometry: Single geometry
            - list: List of Shapely geometries
        resolution (int): S2 resolution level [0..28] (0=coarsest, 28=finest)
        predicate (str, optional): Spatial predicate for polygons:
            - 'intersect': Include cells that intersect the polygon
            - 'within': Include cells completely within the polygon
            - 'centroid_within': Include cells whose centroids are within the polygon
            - 'largest_overlap': Include cells with largest overlap with the polygon
        compact (bool, optional): Enable S2 compact mode to reduce cell count for polygons and lines
        topology (bool, optional): Enable topology preserving mode (currently not implemented)
        output_format (str, optional): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str, optional): Output file path. If None, uses default naming
        include_properties (bool, optional): Whether to include original feature properties in output
        **kwargs: Additional arguments passed to geopandas read functions

    Returns:
        dict or str: Output in the specified format. Returns file path if output_path is specified,
        otherwise returns the data directly.

    Raises:
        ValueError: If input data type is not supported or conversion fails
        TypeError: If resolution is not an integer

    Example:
        >>> # Convert from file
        >>> result = vector2s2("cities.shp", 10, predicate="intersect")
        
        >>> # Convert from GeoDataFrame
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file("cities.shp")
        >>> result = vector2s2(gdf, 10, output_format="geojson")
        
        >>> # Convert from GeoJSON dict
        >>> geojson = {"type": "FeatureCollection", "features": [...]}
        >>> result = vector2s2(geojson, 10, compact=True)
    """
    # Process input data directly
    if hasattr(data, "geometry") and hasattr(data, "columns"):
        # GeoDataFrame
        result = geodataframe2s2(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        # DataFrame with geometry column
        result = dataframe2s2(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        # Shapely geometry objects
        result = geometry2s2(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        # GeoJSON dict
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2s2(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        # File path or URL
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2s2(
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
    Convert GeoJSON FeatureCollection to various output formats.

    Args:
        result (dict): GeoJSON FeatureCollection dictionary
        output_format (str): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str, optional): Output file path. If None, uses default naming

    Returns:
        dict or str: Output in the specified format or file path

    Raises:
        ValueError: If output format is not supported

    Example:
        >>> result = {"type": "FeatureCollection", "features": [...]}
        >>> output = convert_to_output_format(result, "geojson", "output.geojson")
        >>> output
        'output.geojson'
    """
    gdf = gpd.GeoDataFrame.from_features(result["features"])
    gdf.set_crs(epsg=4326, inplace=True)
    if output_format.lower() == "geojson":
        if output_path:
            with open(output_path, "w") as f:
                json.dump(result, f, indent=2)
            return output_path
        else:
            return result
    elif output_format.lower() == "gpkg":
        if output_path:
            gdf.to_file(output_path, driver="GPKG")
            return output_path
        else:
            gdf.to_file("vector2s2.gpkg", driver="GPKG")
            return "vector2s2.gpkg"
    elif output_format.lower() == "parquet":
        if output_path:
            gdf.to_parquet(output_path, index=False)
            return output_path
        else:
            gdf.to_parquet("vector2s2.parquet", index=False)
            return "vector2s2.parquet"
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
            gdf.to_file("vector2s2.shp", driver="ESRI Shapefile")
            return "vector2s2.shp"
    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


# --- CLI ---
def vector2s2_cli():
    """
    Command-line interface for vector2s2 conversion.

    This function provides a command-line interface for converting vector data to S2 grid cells.
    It parses command-line arguments and calls the main vector2s2 function.

    Usage:
        python vector2s2.py -i input.shp -r 10 -p intersect -c -f geojson -o output.geojson

    Arguments:
        -i, --input: Input file path or URL
        -r, --resolution: S2 resolution [0..28]
        -p, --predicate: Spatial predicate (intersect, within, centroid_within, largest_overlap)
        -c, --compact: Enable S2 compact mode
        -t, --topology: Enable topology preserving mode
        -np, --no-props: Do not include original feature properties
        -f, --format: Output format (geojson, gpkg, parquet, csv, shapefile)
        -o, --output: Output file path

    Example:
        >>> # Convert shapefile to S2 cells at resolution 10
        >>> # python vector2s2.py -i cities.shp -r 10 -p intersect -f geojson
    """
    parser = argparse.ArgumentParser(description="Convert vector data to S2 grid cells")
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(29),
        metavar="[0-28]",
        help="S2 resolution [0..28] (0=coarsest, 28=finest)",
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
        help="Enable S2 compact mode for polygons",
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
            args.resolution = validate_s2_resolution(args.resolution)
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
        output_path = f"vector2s2{extensions.get(args.format, '')}"
    try:
        vector2s2(
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
    vector2s2_cli()

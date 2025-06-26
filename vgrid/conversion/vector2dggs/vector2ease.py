import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Polygon, box
import pandas as pd
import geopandas as gpd
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.easedggs.constants import levels_specs, geo_crs, ease_crs
from vgrid.utils.easedggs.dggs.grid_addressing import (
    grid_ids_to_geos,
    geos_to_grid_ids,
    geo_polygon_to_grid_ids,
)
from vgrid.conversion.dggscompact import ease_compact
from vgrid.utils.geometry import check_predicate


def validate_ease_resolution(resolution):
    """
    Validate that EASE resolution is in the valid range [0..6].

    Args:
        resolution (int): Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [0..6]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 0 or resolution > 6:
        raise ValueError(f"Resolution must be in range [0..6], got {resolution}")

    return resolution


def point2ease(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a single point geometry to EASE grid cell.

    Args:
        resolution (int): EASE resolution level [0..6]
        point (shapely.geometry.Point): Point geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for points)
        compact (bool, optional): Enable EASE compact mode (not used for points)
        topology (bool, optional): Enable topology preserving mode (not used for points)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing EASE cells containing the point
    """
    ease_features = []
    latitude = point.y
    longitude = point.x
    ease_cell = geos_to_grid_ids([(longitude, latitude)], level=resolution)
    ease_id = ease_cell["result"]["data"][0]
    level = int(ease_id[1])
    level_spec = levels_specs[level]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]
    geo = grid_ids_to_geos([ease_id])
    center_lon, center_lat = geo["result"]["data"][0]
    cell_min_lat = center_lat - (180 / (2 * n_row))
    cell_max_lat = center_lat + (180 / (2 * n_row))
    cell_min_lon = center_lon - (360 / (2 * n_col))
    cell_max_lon = center_lon + (360 / (2 * n_col))
    cell_polygon = Polygon(
        [
            [cell_min_lon, cell_min_lat],
            [cell_max_lon, cell_min_lat],
            [cell_max_lon, cell_max_lat],
            [cell_min_lon, cell_max_lat],
            [cell_min_lon, cell_min_lat],
        ]
    )
    if cell_polygon:
        num_edges = 4
        ease_feature = geodesic_dggs_to_feature(
            "ease", ease_id, level, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            ease_feature["properties"].update(feature_properties)
        ease_features.append(ease_feature)
    return ease_features


def polyline2ease(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert line geometries (LineString, MultiLineString) to EASE grid cells.

    Args:
        resolution (int): EASE resolution level [0..6]
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Line geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for lines)
        compact (bool, optional): Enable EASE compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (not used for lines)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing EASE cells intersecting the line
    """
    ease_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        poly_bbox = box(*polyline.bounds)
        polygon_bbox_wkt = poly_bbox.wkt
        cells_bbox = geo_polygon_to_grid_ids(
            polygon_bbox_wkt,
            resolution,
            geo_crs,
            ease_crs,
            levels_specs,
            return_centroids=True,
            wkt_geom=True,
        )
        ease_cells = cells_bbox["result"]["data"]
        if compact:
            ease_cells = ease_compact(ease_cells)
        for ease_cell in ease_cells:
            cell_resolution = int(ease_cell[1])
            level_spec = levels_specs[cell_resolution]
            n_row = level_spec["n_row"]
            n_col = level_spec["n_col"]
            geo = grid_ids_to_geos([ease_cell])
            center_lon, center_lat = geo["result"]["data"][0]
            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))
            cell_polygon = Polygon(
                [
                    [cell_min_lon, cell_min_lat],
                    [cell_max_lon, cell_min_lat],
                    [cell_max_lon, cell_max_lat],
                    [cell_min_lon, cell_max_lat],
                    [cell_min_lon, cell_min_lat],
                ]
            )
            if cell_polygon.intersects(polyline):
                num_edges = 4
                ease_feature = geodesic_dggs_to_feature(
                    "ease", str(ease_cell), cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    ease_feature["properties"].update(feature_properties)
                ease_features.append(ease_feature)
    return ease_features


def polygon2ease(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert polygon geometries (Polygon, MultiPolygon) to EASE grid cells.

    Args:
        resolution (int): EASE resolution level [0..6]
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable EASE compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (not used for polygons)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing EASE cells based on predicate
    """
    ease_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        poly_bbox = box(*polygon.bounds)
        polygon_bbox_wkt = poly_bbox.wkt
        cells_bbox = geo_polygon_to_grid_ids(
            polygon_bbox_wkt,
            resolution,
            geo_crs,
            ease_crs,
            levels_specs,
            return_centroids=True,
            wkt_geom=True,
        )
        ease_cells = cells_bbox["result"]["data"]
        if compact:
            ease_cells = ease_compact(ease_cells)
        for ease_cell in ease_cells:
            cell_resolution = int(ease_cell[1])
            level_spec = levels_specs[cell_resolution]
            n_row = level_spec["n_row"]
            n_col = level_spec["n_col"]
            geo = grid_ids_to_geos([ease_cell])
            center_lon, center_lat = geo["result"]["data"][0]
            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))
            cell_polygon = Polygon(
                [
                    [cell_min_lon, cell_min_lat],
                    [cell_max_lon, cell_min_lat],
                    [cell_max_lon, cell_max_lat],
                    [cell_min_lon, cell_max_lat],
                    [cell_min_lon, cell_min_lat],
                ]
            )
            if check_predicate(cell_polygon, polygon, predicate):
                num_edges = 4
                ease_feature = geodesic_dggs_to_feature(
                    "ease", str(ease_cell), cell_resolution, cell_polygon, num_edges
                )
                if feature_properties:
                    ease_feature["properties"].update(feature_properties)
                ease_features.append(ease_feature)
    return ease_features


# --- Main geometry conversion ---
def geometry2ease(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a list of geometries to EASE grid cells.

    Args:
        geometries (shapely.geometry.BaseGeometry or list): Single geometry or list of geometries
        resolution (int): EASE resolution level [0..6]
        properties_list (list, optional): List of property dictionaries for each geometry
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable EASE compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with EASE grid cells
    """
    resolution = validate_ease_resolution(resolution)
    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    ease_features = []
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
            ease_features.extend(point2ease(resolution, geom, props))
        elif geom.geom_type == "MultiPoint":
            for pt in geom.geoms:
                ease_features.extend(point2ease(resolution, pt, props))
        elif geom.geom_type in ("LineString", "MultiLineString"):
            ease_features.extend(
                polyline2ease(
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
            ease_features.extend(
                polygon2ease(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
    return {"type": "FeatureCollection", "features": ease_features}


# --- DataFrame/GeoDataFrame conversion ---
def dataframe2ease(
    df,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a pandas DataFrame with geometry column to EASE grid cells.

    Args:
        df (pandas.DataFrame): DataFrame with geometry column
        resolution (int): EASE resolution level [0..6]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable EASE compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with EASE grid cells
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
    return geometry2ease(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2ease(
    gdf,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a GeoDataFrame to EASE grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int): EASE resolution level [0..6]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable EASE compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with EASE grid cells
    """
    geometries = []
    properties_list = []
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is not None:
            geometries.append(geom)
            props = row.to_dict()
            if "geometry" in props:
                del props["geometry"]
            properties_list.append(props)
    return geometry2ease(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


# --- Main vector2ease function ---
def vector2ease(
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
    Convert vector data to EASE grid cells from various input formats.

    Args:
        data: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
        resolution (int): EASE resolution [0..6]
        compact (bool): Enable EASE compact mode for polygons (default: False)
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
        result = geodataframe2ease(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        # DataFrame with geometry column
        result = dataframe2ease(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        # Shapely geometry objects
        result = geometry2ease(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        # GeoJSON dict
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2ease(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        # File path or URL
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2ease(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")
    return convert_to_output_format(result, output_format, output_path)


# --- Output format conversion ---
def convert_to_output_format(
    result, output_format, output_path=None
):  # topology=Falset to bypass warning from pyogrio
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
            gdf.to_file("vector2ease.gpkg", driver="GPKG")
            return "vector2ease.gpkg"
    elif output_format.lower() == "parquet":
        if output_path:
            gdf.to_parquet(output_path, index=False)
            return output_path
        else:
            return gdf.to_parquet(index=False)
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
            gdf.to_file("vector2ease.shp", driver="ESRI Shapefile")
            return "vector2ease.shp"
    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


# --- CLI ---
def vector2ease_cli():
    """
    Command-line interface for vector2ease conversion.

    This function provides a command-line interface for converting vector data to EASE grid cells.
    It parses command-line arguments and calls the main vector2ease function.

    Usage:
        python vector2ease.py -i input.shp -r 3 -c -f geojson -o output.geojson

    Arguments:
        -i, --input: Input file path or URL
        -r, --resolution: EASE resolution [0..6]
        -c, --compact: Enable EASE compact mode
        -p, --predicate: Spatial predicate (intersect, within, centroid_within, largest_overlap)
        -t, --topology: Enable topology preserving mode
        -np, --no-props: Do not include original feature properties
        -f, --format: Output format (geojson, gpkg, parquet, csv, shapefile)
        -o, --output: Output file path
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to EASE grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(7),
        metavar="[0-6]",
        help="EASE resolution [0..6] (0=coarsest, 6=finest)",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable EASE compact mode for polygons",
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
            args.resolution = validate_ease_resolution(args.resolution)
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
        output_path = f"vector2ease{extensions.get(args.format, '')}"
    try:
        vector2ease(
            data,
            args.resolution,
            predicate=args.predicate,
            topology=args.topology,
            compact=args.compact,
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
    vector2ease_cli()

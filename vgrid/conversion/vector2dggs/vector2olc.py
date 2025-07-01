import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Polygon, shape, MultiPoint, MultiLineString, MultiPolygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import olc
from vgrid.generator.olcgrid import generate_grid, refine_cell
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.dggscompact import olccompact
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
    shortest_polyline_distance,
    shortest_polygon_distance,
)
from vgrid.stats.olcstats import olc_metrics
from math import sqrt

def validate_olc_resolution(resolution):
    """
    Validate that OLC resolution is in the valid range [2,4,6,8,10,11,12,13,14,15].

    Args:
        resolution (int): Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [2,4,6,8,10,11,12,13,14,15]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution not in [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]:
        raise ValueError(
            f"Resolution must be in [2,4,6,8,10,11,12,13,14,15], got {resolution}"
        )

    return resolution


def point2olc(
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
    Convert a point geometry to an OLC grid cell.

    Args:
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        point (shapely.geometry.Point): Point geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for points)
        compact (bool, optional): Enable OLC compact mode (not used for points)
        topology (bool, optional): Enable topology preserving mode - ensures disjoint points have disjoint OLC cells
        include_properties (bool, optional): Whether to include properties in output
        all_points: List of all points for topology preservation (required when topology=True)

    Returns:
        list: List of GeoJSON feature dictionaries representing OLC cells containing the point
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_points is None:
            raise ValueError("all_points parameter is required when topology=True")
        
        # Calculate the shortest distance between all points
        shortest_distance = shortest_point_distance(all_points)
        
        # Find resolution where OLC cell size is smaller than shortest distance
        # This ensures disjoint points have disjoint OLC cells
        if shortest_distance > 0:
            for res in [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]:  # OLC valid resolutions
                _, avg_edge_length, _ = olc_metrics(res)
                # Use a factor to ensure sufficient separation (cell diameter is ~2x edge length)
                cell_diameter = avg_edge_length * sqrt(2)* 2
                if cell_diameter < shortest_distance:
                    resolution = res
                    break
            else:
                # If no resolution found, use the highest resolution
                resolution = 15
        else:
            # Single point or no distance, use provided resolution
            pass

    olc_features = []
    olc_id = olc.encode(point.y, point.x, resolution)
    coord = olc.decode(olc_id)
    if coord:
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi
        cell_polygon = Polygon(
            [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat],
            ]
        )
        olc_feature = graticule_dggs_to_feature("olc", olc_id, resolution, cell_polygon)
        if include_properties and feature_properties:
            olc_feature["properties"].update(feature_properties)
        olc_features.append(olc_feature)
    return olc_features


def polyline2olc(
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
    Convert line geometries (LineString, MultiLineString) to OLC grid cells.

    Args:
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Line geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for lines)
        compact (bool, optional): Enable OLC compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode - ensures disjoint polylines have disjoint OLC cells
        include_properties (bool, optional): Whether to include properties in output
        all_polylines: List of all polylines for topology preservation (required when topology=True)

    Returns:
        list: List of GeoJSON feature dictionaries representing OLC cells intersecting the line
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_polylines is None:
            raise ValueError("all_polylines parameter is required when topology=True")
        
        # Calculate the shortest distance between all polylines
        shortest_distance = shortest_polyline_distance(all_polylines)
        
        # Find resolution where OLC cell size is smaller than shortest distance
        # This ensures disjoint polylines have disjoint OLC cells
        if shortest_distance > 0:
            for res in [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]:  # OLC valid resolutions
                _, avg_edge_length, _ = olc_metrics(res)
                # Use a factor to ensure sufficient separation (cell diameter is ~2x edge length)
                cell_diameter = avg_edge_length * sqrt(2) * 4  # in case there are 2 cells representing the same line segment
                if cell_diameter < shortest_distance:
                    resolution = res
                    break
            else:
                # If no resolution found, use the highest resolution
                resolution = 15
        else:
            # Single polyline or no distance, use provided resolution
            pass
    olc_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        base_resolution = 2
        base_cells = generate_grid(base_resolution, verbose=False)
        seed_cells = []
        for base_cell in base_cells["features"]:
            base_cell_poly = Polygon(base_cell["geometry"]["coordinates"][0])
            if polyline.intersects(base_cell_poly):
                seed_cells.append(base_cell)
        refined_features = []
        for seed_cell in seed_cells:
            seed_cell_poly = Polygon(seed_cell["geometry"]["coordinates"][0])
            if seed_cell_poly.contains(polyline) and resolution == base_resolution:
                refined_features.append(seed_cell)
            else:
                refined_features.extend(
                    refine_cell(
                        seed_cell_poly.bounds, base_resolution, resolution, polyline
                    )
                )
        resolution_features = [
            refined_feature
            for refined_feature in refined_features
            if refined_feature["properties"]["resolution"] == resolution
        ]
        seen_olc_codes = set()
        for resolution_feature in resolution_features:
            olc_id = resolution_feature["properties"]["olc"]
            if olc_id not in seen_olc_codes:
                if include_properties and feature_properties:
                    resolution_feature["properties"].update(feature_properties)
                olc_features.append(resolution_feature)
                seen_olc_codes.add(olc_id)
    olc_geojson = {"type": "FeatureCollection", "features": olc_features}
    if compact:
        return olccompact(olc_geojson)["features"]
    return olc_features


def polygon2olc(
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
    Convert polygon geometries (Polygon, MultiPolygon) to OLC grid cells.

    Args:
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable OLC compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode - ensures disjoint polygons have disjoint OLC cells
        include_properties (bool, optional): Whether to include properties in output
        all_polygons: List of all polygons for topology preservation (required when topology=True)

    Returns:
        list: List of GeoJSON feature dictionaries representing OLC cells based on predicate
    """
    # If topology preservation is enabled, calculate appropriate resolution
    if topology:
        if all_polygons is None:
            raise ValueError("all_polygons parameter is required when topology=True")
        
        # Calculate the shortest distance between all polygons
        shortest_distance = shortest_polygon_distance(all_polygons)
        
        # Find resolution where OLC cell size is smaller than shortest distance
        # This ensures disjoint polygons have disjoint OLC cells
        if shortest_distance > 0:
            for res in [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]:  # OLC valid resolutions
                _, avg_edge_length, _ = olc_metrics(res)
                # Use a factor to ensure sufficient separation (cell diameter is ~2x edge length)
                cell_diameter = avg_edge_length * sqrt(2) * 4
                if cell_diameter < shortest_distance:
                    resolution = res
                    break
            else:
                # If no resolution found, use the highest resolution
                resolution = 15
        else:
            # Single polygon or no distance, use provided resolution
            pass
    olc_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        base_resolution = 2
        base_cells = generate_grid(base_resolution, verbose=False)
        seed_cells = []
        for base_cell in base_cells["features"]:
            base_cell_poly = Polygon(base_cell["geometry"]["coordinates"][0])
            if polygon.intersects(base_cell_poly):
                seed_cells.append(base_cell)
        refined_features = []
        for seed_cell in seed_cells:
            seed_cell_poly = Polygon(seed_cell["geometry"]["coordinates"][0])
            if seed_cell_poly.contains(polygon) and resolution == base_resolution:
                refined_features.append(seed_cell)
            else:
                refined_features.extend(
                    refine_cell(
                        seed_cell_poly.bounds, base_resolution, resolution, polygon
                    )
                )
        resolution_features = [
            refined_feature
            for refined_feature in refined_features
            if refined_feature["properties"]["resolution"] == resolution
        ]
        seen_olc_codes = set()
        for resolution_feature in resolution_features:
            olc_id = resolution_feature["properties"]["olc"]
            if olc_id not in seen_olc_codes:
                cell_geom = shape(resolution_feature["geometry"])
                if not check_predicate(cell_geom, polygon, predicate):
                    continue
                if include_properties and feature_properties:
                    resolution_feature["properties"].update(feature_properties)
                olc_features.append(resolution_feature)
                seen_olc_codes.add(olc_id)
    olc_geojson = {"type": "FeatureCollection", "features": olc_features}
    if compact:
        return olccompact(olc_geojson)["features"]
    return olc_features


def geometry2olc(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a list of geometries to OLC grid cells.

    Args:
        geometries (shapely.geometry.BaseGeometry or list): Single geometry or list of geometries
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        properties_list (list, optional): List of property dictionaries for each geometry
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable OLC compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with OLC grid cells
    """
    resolution = validate_olc_resolution(resolution)

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

    olc_features = []

    for i, geom in enumerate(tqdm(geometries, desc="Processing features")):
        if geom is None:
            continue

        # Get properties for this geometry
        props = properties_list[i] if i < len(properties_list) else {}

        # Process based on geometry type
        if geom.geom_type == "Point":
            point_features = point2olc(
                resolution,
                geom,
                props,
                predicate,
                compact,
                topology,
                include_properties,
                all_points,  # Pass all points for topology preservation
            )
            olc_features.extend(point_features)

        elif geom.geom_type == "MultiPoint":
            for point in geom.geoms:
                point_features = point2olc(
                    resolution,
                    point,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                    all_points,  # Pass all points for topology preservation
                )
                olc_features.extend(point_features)

        elif geom.geom_type in ["LineString", "MultiLineString"]:
            polyline_features = polyline2olc(
                resolution,
                geom,
                props,
                predicate,
                compact,
                topology,
                include_properties,
                all_polylines,  # Pass all polylines for topology preservation
            )
            olc_features.extend(polyline_features)

        elif geom.geom_type in ["Polygon", "MultiPolygon"]:
            poly_features = polygon2olc(
                resolution,
                geom,
                props,
                predicate,
                compact,
                topology,
                include_properties,
                all_polygons=all_polygons,  # Pass all polygons for topology preservation
            )
            olc_features.extend(poly_features)

        else:
            raise ValueError(f"Unsupported geometry type: {geom.geom_type}")

    return {"type": "FeatureCollection", "features": olc_features}


def dataframe2olc(
    df,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a pandas DataFrame with geometry column to OLC grid cells.

    Args:
        df (pandas.DataFrame): DataFrame with geometry column
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable OLC compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with OLC grid cells
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
    return geometry2olc(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2olc(
    gdf,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a GeoDataFrame to OLC grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        predicate (str, optional): Spatial predicate to apply for polygons
        compact (bool, optional): Enable OLC compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        dict: GeoJSON FeatureCollection with OLC grid cells
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
    return geometry2olc(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def vector2olc(
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
    Convert vector data to OLC grid cells from various input formats.

    Args:
        data: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
        resolution (int): OLC resolution [2,4,6,8,10,11,12,13,14,15]
        compact (bool): Enable OLC compact mode for polygons (default: False)
        output_format (str): Output format ('geojson', 'gpkg', 'parquet', 'csv', 'shapefile')
        output_path (str): Output file path (optional)
        include_properties (bool): If False, do not include original feature properties. (default: True)
        **kwargs: Additional arguments passed to geopandas read functions
    Returns:
        dict or str: Output in the specified format
    """
    if hasattr(data, "geometry") and hasattr(data, "columns"):
        result = geodataframe2olc(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        result = dataframe2olc(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        result = geometry2olc(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2olc(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2olc(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")
    return convert_to_output_format(result, output_format, output_path)


def convert_to_output_format(result, output_format, output_path=None):
    """
    Convert OLC result to specified output format.

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

        # Set CRS to WGS84 (EPSG:4326) since OLC uses WGS84 coordinates
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
            gdf.to_file("vector2olc.gpkg", driver="GPKG")
            return "vector2olc.gpkg"

    elif output_format.lower() == "parquet":
        if output_path:
            gdf.to_parquet(output_path, index=False)
            return output_path
        else:
            gdf.to_parquet("vector2olc.parquet", index=False)
            return "vector2olc.parquet"

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
            gdf.to_file("vector2olc.shp", driver="ESRI Shapefile")
            return "vector2olc.shp"

    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


def vector2olc_cli():
    """
    Command-line interface for vector2olc conversion.

    Usage:
        python vector2olc.py -i input.shp -r 10 -c -f geojson -o output.geojson

    Arguments:
        -i, --input: Input file path or URL
        -r, --resolution: OLC resolution (see OLC spec)
        -c, --compact: Enable OLC compact mode
        -p, --predicate: Spatial predicate (intersect, within, centroid_within, largest_overlap)
        -t, --topology: Enable topology preserving mode
        -np, --no-props: Do not include original feature properties
        -f, --format: Output format (geojson, gpkg, parquet, csv, shapefile)
        -o, --output: Output file path
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to OLC grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=[2, 4, 6, 8, 10, 11, 12, 13, 14, 15],
        metavar="[2,4,6,8,10,11,12,13,14,15]",
        help="OLC resolution (see OLC spec)",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable OLC compact mode for polygons",
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate: intersect, within, centroid_within, largest_overlap for polygons",
    )
    parser.add_argument(
        "-t", "--topology", action="store_true", 
        help="Enable topology preserving mode ensuring disjoint features have disjoint OLC cells by automatically calculating appropriate resolution."
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
            args.resolution = validate_olc_resolution(args.resolution)
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
        output_path = f"vector2olc{extensions.get(args.format, '')}"
    try:
        vector2olc(
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
    vector2olc_cli()

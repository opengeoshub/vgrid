import sys
import argparse
import json
import platform
from tqdm import tqdm
from shapely.geometry import box
import pandas as pd
import geopandas as gpd
from vgrid.utils.geometry import check_predicate


def is_windows():
    return platform.system() == "Windows"


if is_windows():
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.generator.isea3hgrid import (
        isea3h_cell_to_polygon,
        isea3h_res_accuracy_dict,
        isea3h_accuracy_res_dict,
        get_isea3h_children_cells_within_bbox,
    )
    from vgrid.conversion.dggscompact import isea3h_compact
    from vgrid.generator.settings import geodesic_dggs_to_feature


def validate_isea3h_resolution(resolution):
    """
    Validate that ISEA3H resolution is in the valid range [0..32].

    Args:
        resolution: Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [0..32]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

        if resolution < 0 or resolution > 32:
            raise ValueError(f"Resolution must be in range [0..32], got {resolution}")

    return resolution


def point2isea3h(
    isea3h_dggs,
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    isea3h_features = []
    accuracy = isea3h_res_accuracy_dict.get(resolution)
    lat_long_point = LatLongPoint(point.y, point.x, accuracy)
    isea3h_cell = isea3h_dggs.convert_point_to_dggs_cell(lat_long_point)
    cell_polygon = isea3h_cell_to_polygon(isea3h_dggs, isea3h_cell)
    if cell_polygon:
        isea3h_id = isea3h_cell.get_cell_id()
        cell_resolution = resolution
        num_edges = 3 if cell_resolution == 0 else 6
        isea3h_feature = geodesic_dggs_to_feature(
            "isea3h", isea3h_id, cell_resolution, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            isea3h_feature["properties"].update(feature_properties)
        isea3h_features.append(isea3h_feature)
    return isea3h_features


def polyline2isea3h(
    isea3h_dggs,
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    isea3h_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        accuracy = isea3h_res_accuracy_dict.get(resolution)
        bounding_box = box(*polyline.bounds)
        bounding_box_wkt = bounding_box.wkt
        shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(
            bounding_box_wkt, ShapeStringFormat.WKT, accuracy
        )
        shape = shapes[0]
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_child_cells = get_isea3h_children_cells_within_bbox(
            isea3h_dggs, bounding_cell.get_cell_id(), bounding_box, resolution
        )
        if compact:
            bounding_child_cells = isea3h_compact(isea3h_dggs, bounding_child_cells)
        for child in bounding_child_cells:
            isea3h_cell = DggsCell(child)
            cell_polygon = isea3h_cell_to_polygon(isea3h_dggs, isea3h_cell)
            if cell_polygon.intersects(polyline):
                isea3h_id = isea3h_cell.get_cell_id()
                isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)
                cell_accuracy = isea3h2point._accuracy
                cell_resolution = isea3h_accuracy_res_dict.get(cell_accuracy)
                num_edges = 3 if cell_resolution == 0 else 6
                isea3h_feature = geodesic_dggs_to_feature(
                    "isea3h", isea3h_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    isea3h_feature["properties"].update(feature_properties)
                isea3h_features.append(isea3h_feature)
    return isea3h_features


def polygon2isea3h(
    isea3h_dggs,
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    isea3h_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        accuracy = isea3h_res_accuracy_dict.get(resolution)
        bounding_box = box(*polygon.bounds)
        bounding_box_wkt = bounding_box.wkt
        shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(
            bounding_box_wkt, ShapeStringFormat.WKT, accuracy
        )
        shape = shapes[0]
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_child_cells = get_isea3h_children_cells_within_bbox(
            isea3h_dggs, bounding_cell.get_cell_id(), bounding_box, resolution
        )
        if compact:
            bounding_child_cells = isea3h_compact(isea3h_dggs, bounding_child_cells)
        for child in bounding_child_cells:
            isea3h_cell = DggsCell(child)
            cell_polygon = isea3h_cell_to_polygon(isea3h_dggs, isea3h_cell)
            if check_predicate(cell_polygon, polygon, predicate):
                isea3h_id = isea3h_cell.get_cell_id()
                isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)
                cell_accuracy = isea3h2point._accuracy
                cell_resolution = isea3h_accuracy_res_dict.get(cell_accuracy)
                num_edges = 3 if cell_resolution == 0 else 6
                isea3h_feature = geodesic_dggs_to_feature(
                    "isea3h", isea3h_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    isea3h_feature["properties"].update(feature_properties)
                isea3h_features.append(isea3h_feature)
    return isea3h_features


def geometry2isea3h(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    if not is_windows():
        raise NotImplementedError("ISEA3H DGGS conversion is only supported on Windows")

    resolution = validate_isea3h_resolution(resolution)
    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    isea3h_dggs = Eaggr(Model.ISEA3H)
    isea3h_features = []
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
            isea3h_features.extend(
                point2isea3h(
                    isea3h_dggs,
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
                isea3h_features.extend(point2isea3h(isea3h_dggs, resolution, pt, props))
        elif geom.geom_type in ("LineString", "MultiLineString"):
            isea3h_features.extend(
                polyline2isea3h(
                    isea3h_dggs,
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
            isea3h_features.extend(
                polygon2isea3h(
                    isea3h_dggs,
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
    return {"type": "FeatureCollection", "features": isea3h_features}


def dataframe2isea3h(
    df,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
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
    return geometry2isea3h(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2isea3h(
    gdf,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
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
    return geometry2isea3h(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def vector2isea3h(
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
    if not is_windows():
        raise NotImplementedError("ISEA3H DGGS conversion is only supported on Windows")

    if hasattr(data, "geometry") and hasattr(data, "columns"):
        result = geodataframe2isea3h(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        result = dataframe2isea3h(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        result = geometry2isea3h(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2isea3h(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2isea3h(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to read file/URL {data}: {str(e)}")
    else:
        raise ValueError(f"Unsupported input type: {type(data)}")
    return convert_to_output_format(result, output_format, output_path)


def convert_to_output_format(result, output_format, output_path=None):
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
            gdf.to_file("vector2isea3h.gpkg", driver="GPKG")
            return "vector2isea3h.gpkg"
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
            gdf.to_file("vector2isea3h.shp", driver="ESRI Shapefile")
            return "vector2isea3h.shp"
    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


def vector2isea3h_cli():
    parser = argparse.ArgumentParser(
        description="Convert vector data to ISEA3H grid cells (Windows only)"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(33),
        metavar="[0-32]",
        help="ISEA3H resolution [0..32] (0=coarsest, 32=finest)",
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
        "-c",
        "--compact",
        action="store_true",
        help="Enable ISEA3H compact mode for polygons",
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
    if not is_windows():
        print("ISEA3H DGGS conversion is only supported on Windows", file=sys.stderr)
        sys.exit(1)
    # Validate resolution if provided
    if args.resolution is not None:
        try:
            args.resolution = validate_isea3h_resolution(args.resolution)
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
        output_path = f"vector2isea3h{extensions.get(args.format, '')}"
    try:
        vector2isea3h(
            data,
            args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            include_properties=args.include_properties,
            output_format=args.format,
            output_path=output_path,
        )
        if output_path:
            print(f"Output saved to {output_path}")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2isea3h_cli()

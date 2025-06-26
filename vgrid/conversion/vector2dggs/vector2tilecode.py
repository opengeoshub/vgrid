import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Polygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import tilecode
from vgrid.utils import mercantile
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.dggscompact import tilecodecompact
import re
from vgrid.utils.geometry import check_predicate


def validate_tilecode_resolution(resolution):
    """
    Validate that Tilecode resolution is in the valid range [0..29] (0=coarsest, 29=finest).

    Args:
        resolution: Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [0..29]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 0 or resolution > 29:
        raise ValueError(f"Resolution must be in range [0..29], got {resolution}")

    return resolution


def point2tilecode(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    tilecode_features = []
    tilecode_id = tilecode.latlon2tilecode(point.y, point.x, resolution)
    tilecode_cell = mercantile.tile(point.x, point.y, resolution)
    bounds = mercantile.bounds(tilecode_cell)
    if bounds:
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east
        cell_polygon = Polygon(
            [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat],
            ]
        )
        tilecode_feature = graticule_dggs_to_feature(
            "tilecode", tilecode_id, resolution, cell_polygon
        )
        if include_properties and feature_properties:
            tilecode_feature["properties"].update(feature_properties)
        tilecode_features.append(tilecode_feature)
    return tilecode_features


def polyline2tilecode(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    tilecode_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        tilecodes = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        tilecode_ids = []
        for tile in tilecodes:
            tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
            tilecode_ids.append(tilecode_id)
        for tilecode_id in tilecode_ids:
            match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
            if not match:
                raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")
            cell_resolution = int(match.group(1))
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                cell_polygon = Polygon(
                    [
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat],
                    ]
                )
                if cell_polygon.intersects(polyline):
                    tilecode_feature = graticule_dggs_to_feature(
                        "tilecode", tilecode_id, cell_resolution, cell_polygon
                    )
                    if feature_properties:
                        tilecode_feature["properties"].update(feature_properties)
                    tilecode_features.append(tilecode_feature)
    if compact:
        tilecode_geojson = {"type": "FeatureCollection", "features": tilecode_features}
        return tilecodecompact(tilecode_geojson)["features"]
    return tilecode_features


def polygon2tilecode(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    tilecode_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        tilecodes = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        tilecode_ids = []
        for tile in tilecodes:
            tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
            tilecode_ids.append(tilecode_id)
        for tilecode_id in tilecode_ids:
            match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
            if not match:
                raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")
            cell_resolution = int(match.group(1))
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                cell_polygon = Polygon(
                    [
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat],
                    ]
                )
                if check_predicate(cell_polygon, polygon, predicate):
                    tilecode_feature = graticule_dggs_to_feature(
                        "tilecode", tilecode_id, cell_resolution, cell_polygon
                    )
                    if include_properties and feature_properties:
                        tilecode_feature["properties"].update(feature_properties)
                    tilecode_features.append(tilecode_feature)
    if compact:
        tilecode_geojson = {"type": "FeatureCollection", "features": tilecode_features}
        return tilecodecompact(tilecode_geojson)["features"]
    return tilecode_features


def geometry2tilecode(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    resolution = validate_tilecode_resolution(resolution)

    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    tilecode_features = []
    for idx, geom in tqdm(enumerate(geometries), desc="Processing features"):
        props = (
            properties_list[idx]
            if properties_list and idx < len(properties_list)
            else {}
        )
        if not include_properties:
            props = {}
        if geom.geom_type == "Point":
            tilecode_features.extend(
                point2tilecode(
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
                tilecode_features.extend(
                    point2tilecode(
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
            tilecode_features.extend(
                polyline2tilecode(
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
            tilecode_features.extend(
                polygon2tilecode(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
    return {"type": "FeatureCollection", "features": tilecode_features}


def dataframe2tilecode(
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
    return geometry2tilecode(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2tilecode(
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
    return geometry2tilecode(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def vector2tilecode(
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
    if hasattr(data, "geometry") and hasattr(data, "columns"):
        result = geodataframe2tilecode(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        result = dataframe2tilecode(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        result = geometry2tilecode(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2tilecode(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2tilecode(
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
            gdf.to_file("vector2tilecode.gpkg", driver="GPKG")
            return "vector2tilecode.gpkg"
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
            gdf.to_file("vector2tilecode.shp", driver="ESRI Shapefile")
            return "vector2tilecode.shp"
    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


def vector2tilecode_cli():
    parser = argparse.ArgumentParser(
        description="Convert vector data to Tilecode grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(0, 30),
        metavar="[0-29]",
        help="Tilecode resolution [0..29] (0=coarsest, 29=finest)",
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
        help="Enable Tilecode compact mode for polygons",
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
            args.resolution = validate_tilecode_resolution(args.resolution)
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
        output_path = f"vector2tilecode{extensions.get(args.format, '')}"
    try:
        vector2tilecode(
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
    vector2tilecode_cli()

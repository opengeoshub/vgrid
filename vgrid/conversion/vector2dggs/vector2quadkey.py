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
from vgrid.conversion.dggscompact import quadkeycompact
from vgrid.utils.geometry import check_predicate


def validate_quadkey_resolution(resolution):
    """
    Validate that Quadkey resolution is in the valid range [0..29] (0=coarsest, 29=finest).

    Args:
        resolution: Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [0..29] (0=coarsest, 29=finest)
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 0 or resolution > 29:
        raise ValueError(
            f"Resolution must be in range [0..29] (0=coarsest, 29=finest), got {resolution}"
        )

    return resolution


def point2quadkey(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    quadkey_features = []
    quadkey_id = tilecode.latlon2quadkey(point.y, point.x, resolution)
    quadkey_cell = mercantile.tile(point.x, point.y, resolution)
    bounds = mercantile.bounds(quadkey_cell)
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
        quadkey_feature = graticule_dggs_to_feature(
            "quadkey", quadkey_id, resolution, cell_polygon
        )
        if feature_properties:
            quadkey_feature["properties"].update(feature_properties)
        quadkey_features.append(quadkey_feature)
    return quadkey_features


def polyline2quadkey(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    quadkey_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                quadkey_id = mercantile.quadkey(tile)
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
                    quadkey_feature = graticule_dggs_to_feature(
                        "quadkey", quadkey_id, resolution, cell_polygon
                    )
                    if feature_properties:
                        quadkey_feature["properties"].update(feature_properties)
                    quadkey_features.append(quadkey_feature)

    quadkey_geojson = {"type": "FeatureCollection", "features": quadkey_features}

    if compact:
        return quadkeycompact(quadkey_geojson)["features"]
    else:
        return quadkey_features


def polygon2quadkey(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    quadkey_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                quadkey_id = mercantile.quadkey(tile)
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
                    quadkey_feature = graticule_dggs_to_feature(
                        "quadkey", quadkey_id, resolution, cell_polygon
                    )
                    if feature_properties:
                        quadkey_feature["properties"].update(feature_properties)
                    quadkey_features.append(quadkey_feature)

    quadkey_geojson = {"type": "FeatureCollection", "features": quadkey_features}

    if compact:
        return quadkeycompact(quadkey_geojson)["features"]
    else:
        return quadkey_features


def geometry2quadkey(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    resolution = validate_quadkey_resolution(resolution)

    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    quadkey_features = []
    for idx, geom in tqdm(enumerate(geometries), desc="Processing features"):
        props = (
            properties_list[idx]
            if properties_list and idx < len(properties_list)
            else {}
        )
        if not include_properties:
            props = {}
        if geom.geom_type == "Point":
            quadkey_features.extend(
                point2quadkey(
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
                quadkey_features.extend(
                    point2quadkey(
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
            quadkey_features.extend(
                polyline2quadkey(
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
            quadkey_features.extend(
                polygon2quadkey(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
    return {"type": "FeatureCollection", "features": quadkey_features}


def dataframe2quadkey(
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

    return geometry2quadkey(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2quadkey(
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
    return geometry2quadkey(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def vector2quadkey(
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
        result = geodataframe2quadkey(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        result = dataframe2quadkey(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        result = geometry2quadkey(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2quadkey(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2quadkey(
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
            gdf.to_file("vector2quadkey.gpkg", driver="GPKG")
            return "vector2quadkey.gpkg"
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
            gdf.to_file("vector2quadkey.shp", driver="ESRI Shapefile")
            return "vector2quadkey.shp"
    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


def vector2quadkey_cli():
    parser = argparse.ArgumentParser(
        description="Convert vector data to Quadkey grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(0, 30),
        metavar="[0-29]",
        help="Quadkey resolution [0..29] (0=coarsest, 29=finest)",
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
        help="Enable Quadkey compact mode for polygons",
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
            args.resolution = validate_quadkey_resolution(args.resolution)
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
        output_path = f"vector2quadkey{extensions.get(args.format, '')}"
    try:
        vector2quadkey(
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
    vector2quadkey_cli()

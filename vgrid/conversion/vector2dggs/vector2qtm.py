import sys
import argparse
import json
from tqdm import tqdm
from shapely.geometry import Polygon
import pandas as pd
import geopandas as gpd
from vgrid.utils import qtm
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.conversion.dggscompact import qtmcompact
from vgrid.utils.geometry import check_predicate

# QTM facet points
p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (
    (90.0, -180.0),
    (90.0, -90.0),
    (90.0, 0.0),
    (90.0, 90.0),
    (90.0, 180.0),
)
p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (
    (0.0, -180.0),
    (0.0, -90.0),
    (0.0, 0.0),
    (0.0, 90.0),
    (0.0, 180.0),
)
n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (
    (-90.0, -180.0),
    (-90.0, -90.0),
    (-90.0, 0.0),
    (-90.0, 90.0),
    (-90.0, 180.0),
)


def validate_qtm_resolution(resolution):
    """
    Validate that QTM resolution is in the valid range [1..24].

    Args:
        resolution: Resolution value to validate

    Returns:
        int: Validated resolution value

    Raises:
        ValueError: If resolution is not in range [1..24]
        TypeError: If resolution is not an integer
    """
    if not isinstance(resolution, int):
        raise TypeError(
            f"Resolution must be an integer, got {type(resolution).__name__}"
        )

    if resolution < 1 or resolution > 24:
        raise ValueError(f"Resolution must be in range [1..24], got {resolution}")

    return resolution


def point2qtm(
    resolution,
    point,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    qtm_features = []
    latitude = point.y
    longitude = point.x
    qtm_id = qtm.latlon_to_qtm_id(latitude, longitude, resolution)
    facet = qtm.qtm_id_to_facet(qtm_id)
    cell_polygon = qtm.constructGeometry(facet)
    if cell_polygon:
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature(
            "qtm", qtm_id, resolution, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            qtm_feature["properties"].update(feature_properties)
        qtm_features.append(qtm_feature)
    return qtm_features


def polyline2qtm(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    qtm_features = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        levelFacets = {}
        QTMID = {}
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []
            if lvl == 0:
                initial_facets = [
                    [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                    [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                    [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                    [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                    [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                    [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                    [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                    [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
                ]
                for i, facet in enumerate(initial_facets):
                    QTMID[0].append(str(i + 1))
                    levelFacets[0].append(facet)
                    facet_geom = qtm.constructGeometry(facet)
                    if Polygon(facet_geom).intersects(polyline) and resolution == 1:
                        qtm_id = QTMID[0][i]
                        num_edges = 3
                        qtm_feature = geodesic_dggs_to_feature(
                            "qtm", qtm_id, resolution, facet_geom, num_edges
                        )
                        if include_properties and feature_properties:
                            qtm_feature["properties"].update(feature_properties)
                        qtm_features.append(qtm_feature)
                        return qtm_features
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if Polygon(subfacet_geom).intersects(polyline):
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:
                                num_edges = 3
                                qtm_feature = geodesic_dggs_to_feature(
                                    "qtm", new_id, resolution, subfacet_geom, num_edges
                                )
                                if include_properties and feature_properties:
                                    qtm_feature["properties"].update(feature_properties)
                                qtm_features.append(qtm_feature)
    if compact:
        qtm_geojson = {"type": "FeatureCollection", "features": qtm_features}
        return qtmcompact(qtm_geojson)["features"]

    return qtm_features


def polygon2qtm(
    resolution,
    feature,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    qtm_features = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        levelFacets = {}
        QTMID = {}
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []
            if lvl == 0:
                initial_facets = [
                    [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                    [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                    [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                    [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                    [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                    [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                    [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                    [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
                ]
                for i, facet in enumerate(initial_facets):
                    QTMID[0].append(str(i + 1))
                    levelFacets[0].append(facet)
                    facet_geom = qtm.constructGeometry(facet)
                    if Polygon(facet_geom).intersects(polygon) and resolution == 1:
                        qtm_id = QTMID[0][i]
                        num_edges = 3
                        qtm_feature = geodesic_dggs_to_feature(
                            "qtm", qtm_id, resolution, facet_geom, num_edges
                        )
                        if include_properties and feature_properties:
                            qtm_feature["properties"].update(feature_properties)
                        qtm_features.append(qtm_feature)
                        return qtm_features
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if Polygon(subfacet_geom).intersects(polygon):
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:
                                if not check_predicate(
                                    Polygon(subfacet_geom), polygon, predicate
                                ):
                                    continue
                                num_edges = 3
                                qtm_feature = geodesic_dggs_to_feature(
                                    "qtm", new_id, resolution, subfacet_geom, num_edges
                                )
                                if include_properties and feature_properties:
                                    qtm_feature["properties"].update(feature_properties)
                                qtm_features.append(qtm_feature)
    if compact:
        qtm_geojson = {"type": "FeatureCollection", "features": qtm_features}
        return qtmcompact(qtm_geojson)["features"]
    return qtm_features


def geometry2qtm(
    geometries,
    resolution,
    properties_list=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    resolution = validate_qtm_resolution(resolution)

    # Handle single geometry or list of geometries
    if not isinstance(geometries, list):
        geometries = [geometries]

    # Handle properties
    if properties_list is None:
        properties_list = [{} for _ in geometries]
    elif not isinstance(properties_list, list):
        properties_list = [properties_list for _ in geometries]

    qtm_features = []
    for idx, geom in enumerate(tqdm(geometries, desc="Processing features")):
        props = (
            properties_list[idx]
            if properties_list and idx < len(properties_list)
            else {}
        )
        if not include_properties:
            props = {}
        if geom.geom_type == "Point":
            qtm_features.extend(
                point2qtm(
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
                qtm_features.extend(
                    point2qtm(
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
            qtm_features.extend(
                polyline2qtm(
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
            qtm_features.extend(
                polygon2qtm(
                    resolution,
                    geom,
                    props,
                    predicate,
                    compact,
                    topology,
                    include_properties,
                )
            )
    return {"type": "FeatureCollection", "features": qtm_features}


def dataframe2qtm(
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
    return geometry2qtm(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def geodataframe2qtm(
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
    return geometry2qtm(
        geometries,
        resolution,
        properties_list,
        predicate,
        compact,
        topology,
        include_properties,
    )


def vector2qtm(
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
        result = geodataframe2qtm(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif isinstance(data, pd.DataFrame):
        result = dataframe2qtm(
            data, resolution, predicate, compact, topology, include_properties
        )
    elif hasattr(data, "geom_type") or (
        isinstance(data, list) and len(data) > 0 and hasattr(data[0], "geom_type")
    ):
        result = geometry2qtm(
            data, resolution, None, predicate, compact, topology, include_properties
        )
    elif isinstance(data, dict) and "type" in data:
        try:
            gdf = gpd.GeoDataFrame.from_features(data["features"])
            result = geodataframe2qtm(
                gdf, resolution, predicate, compact, topology, include_properties
            )
        except Exception as e:
            raise ValueError(f"Failed to convert GeoJSON to GeoDataFrame: {str(e)}")
    elif isinstance(data, str):
        try:
            gdf = gpd.read_file(data, **kwargs)
            result = geodataframe2qtm(
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
            gdf.to_file("vector2qtm.gpkg", driver="GPKG")
            return "vector2qtm.gpkg"
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
            gdf.to_file("vector2qtm.shp", driver="ESRI Shapefile")
            return "vector2qtm.shp"
    else:
        raise ValueError(
            f"Unsupported output format: {output_format}. Supported formats: geojson, gpkg, parquet, csv, shapefile"
        )


def vector2qtm_cli():
    parser = argparse.ArgumentParser(
        description="Convert vector data to QTM grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(1, 25),
        metavar="[1-24]",
        help="QTM resolution [1..24] (1=coarsest, 24=finest)",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable QTM compact mode for polygons",
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
            args.resolution = validate_qtm_resolution(args.resolution)
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
        output_path = f"vector2qtm{extensions.get(args.format, '')}"
    try:
        vector2qtm(
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
    vector2qtm_cli()

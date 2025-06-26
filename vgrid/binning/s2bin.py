import argparse
import os
import json
import statistics
from collections import defaultdict, Counter
from vgrid.utils import s2
from shapely.geometry import Point, Polygon
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.conversion.latlon2dggs import latlon2s2
from vgrid.utils.antimeridian import fix_polygon


def s2_bin(point_features, resolution, stats, category, field_name):
    s2_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

    for feature in tqdm(point_features, desc="Binning points"):
        geom = feature["geometry"]
        props = feature.get("properties", {})

        if geom["type"] == "Point":
            point = Point(geom["coordinates"])
            s2_token = latlon2s2(point.y, point.x, resolution)
            append_stats_value(s2_bins, s2_token, props, stats, category, field_name)

        elif geom["type"] == "MultiPoint":
            for coords in geom["coordinates"]:
                point = Point(coords)
                s2_token = latlon2s2(point.y, point.x, resolution)
                append_stats_value(
                    s2_bins, s2_token, props, stats, category, field_name
                )

    s2_features = []
    for s2_token, categories in s2_bins.items():
        s2_id = s2.CellId.from_token(s2_token)
        s2_cell = s2.Cell(s2_id)

        vertices = [s2_cell.get_vertex(i) for i in range(4)]
        # Prepare vertices in (longitude, latitude) format for Shapely
        shapely_vertices = []
        for vertex in vertices:
            lat_lng = s2.LatLng.from_point(vertex)  # Convert Point to LatLng
            longitude = lat_lng.lng().degrees  # Access longitude in degrees
            latitude = lat_lng.lat().degrees  # Access latitude in degrees
            shapely_vertices.append((longitude, latitude))

        # Close the polygon by adding the first vertex again
        shapely_vertices.append(shapely_vertices[0])  # Closing the polygon
        # Create a Shapely Polygon
        cell_polygon = fix_polygon(Polygon(shapely_vertices))  # Fix antimeridian

        # polygon = Polygon([(lon, lat) for lat, lon in fix_s2_antimeridian_cells(cell_boundary)])

        if not cell_polygon.is_valid:
            continue
        num_edges = 4
        s2_feature = geodesic_dggs_to_feature(
            "s2", s2_token, resolution, cell_polygon, num_edges
        )

        for cat, values in categories.items():
            key_prefix = "" if category is None else f"{cat}_"

            if stats == "count":
                s2_feature["properties"][f"{key_prefix}count"] = values["count"]
            elif stats == "sum":
                s2_feature["properties"][f"{key_prefix}sum"] = sum(values["sum"])
            elif stats == "min":
                s2_feature["properties"][f"{key_prefix}min"] = min(values["min"])
            elif stats == "max":
                s2_feature["properties"][f"{key_prefix}max"] = max(values["max"])
            elif stats == "mean":
                s2_feature["properties"][f"{key_prefix}mean"] = statistics.mean(
                    values["mean"]
                )
            elif stats == "median":
                s2_feature["properties"][f"{key_prefix}median"] = statistics.median(
                    values["median"]
                )
            elif stats == "std":
                s2_feature["properties"][f"{key_prefix}std"] = (
                    statistics.stdev(values["std"]) if len(values["std"]) > 1 else 0
                )
            elif stats == "var":
                s2_feature["properties"][f"{key_prefix}var"] = (
                    statistics.variance(values["var"]) if len(values["var"]) > 1 else 0
                )
            elif stats == "range":
                s2_feature["properties"][f"{key_prefix}range"] = (
                    max(values["range"]) - min(values["range"])
                    if values["range"]
                    else 0
                )

            elif stats == "minority":
                freq = Counter(values["values"])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                s2_feature["properties"][f"{key_prefix}minority"] = min_item

            elif stats == "majority":
                freq = Counter(values["values"])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                s2_feature["properties"][f"{key_prefix}majority"] = max_item

            elif stats == "variety":
                s2_feature["properties"][f"{key_prefix}variety"] = len(
                    set(values["values"])
                )

        s2_features.append(s2_feature)

    return s2_features


def main():
    parser = argparse.ArgumentParser(description="Binning point to S2 DGGS")
    parser.add_argument(
        "-point",
        "--point",
        type=str,
        required=True,
        help="GeoJSON file path (Point or MultiPoint)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=13,
        help="Resolution of the grid [0..30]",
    )
    parser.add_argument(
        "-stats",
        "--statistics",
        choices=[
            "count",
            "min",
            "max",
            "sum",
            "mean",
            "median",
            "std",
            "var",
            "range",
            "minority",
            "majority",
            "variety",
        ],
        required=True,
        help="Statistic option",
    )

    parser.add_argument(
        "-category",
        "--category",
        required=False,
        help="Optional category field for grouping",
    )
    parser.add_argument(
        "-field", "--field", required=False, help="Numeric field to compute statistics"
    )

    args = parser.parse_args()

    resolution = args.resolution
    point = args.point
    stats = args.statistics
    category = args.category
    field_name = args.field

    if resolution < 0 or resolution > 30:
        print("Error: Please select a resolution in [0..30].")
        return

    if not os.path.exists(point):
        print(f"Error: The file {point} does not exist.")
        return

    if stats != "count" and not field_name:
        print("Error: A field name is required for statistics other than 'count'.")
        return

    with open(point, "r", encoding="utf-8") as f:
        point_data = json.load(f)

    point_features = point_data["features"]
    s2_features = s2_bin(point_features, resolution, stats, category, field_name)

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_bin_s2_{resolution}_{stats}.geojson"

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": s2_features}, f, indent=2)

    print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()

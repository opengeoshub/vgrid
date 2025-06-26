import argparse
import os
import json
import statistics
from collections import defaultdict, Counter
from shapely.geometry import Point
from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.conversion.latlon2dggs import latlon2rhealpix
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.conversion.dggs2geojson import rhealpix_cell_to_polygon


def rhealpix_bin(
    rhealpix_dggs, point_features, resolution, stats, category, field_name
):
    rhealpix_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

    for point_feature in tqdm(point_features, desc="Binning points"):
        geom = point_feature["geometry"]
        props = point_feature.get("properties", {})

        if geom["type"] == "Point":
            point = Point(geom["coordinates"])
            rhealpix_id = latlon2rhealpix(point.y, point.x, resolution)
            append_stats_value(
                rhealpix_bins, rhealpix_id, props, stats, category, field_name
            )

        elif geom["type"] == "MultiPoint":
            for coords in geom["coordinates"]:
                point = Point(coords)
                rhealpix_id = latlon2rhealpix(point.y, point.x, resolution)
                append_stats_value(
                    rhealpix_bins, rhealpix_id, props, stats, category, field_name
                )

    rhealpix_features = []
    for rhealpix_id, categories in rhealpix_bins.items():
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
        num_edges = 4
        if rhealpix_cell.ellipsoidal_shape() == "dart":
            num_edges = 3

        if not cell_polygon.is_valid:
            continue

        rhealpix_feature = geodesic_dggs_to_feature(
            "rhealpix", rhealpix_id, resolution, cell_polygon, num_edges
        )

        for cat, values in categories.items():
            key_prefix = "" if category is None else f"{cat}_"

            if stats == "count":
                rhealpix_feature["properties"][f"{key_prefix}count"] = values["count"]
            elif stats == "sum":
                rhealpix_feature["properties"][f"{key_prefix}sum"] = sum(values["sum"])
            elif stats == "min":
                rhealpix_feature["properties"][f"{key_prefix}min"] = min(values["min"])
            elif stats == "max":
                rhealpix_feature["properties"][f"{key_prefix}max"] = max(values["max"])
            elif stats == "mean":
                rhealpix_feature["properties"][f"{key_prefix}mean"] = statistics.mean(
                    values["mean"]
                )
            elif stats == "median":
                rhealpix_feature["properties"][f"{key_prefix}median"] = (
                    statistics.median(values["median"])
                )
            elif stats == "std":
                rhealpix_feature["properties"][f"{key_prefix}std"] = (
                    statistics.stdev(values["std"]) if len(values["std"]) > 1 else 0
                )
            elif stats == "var":
                rhealpix_feature["properties"][f"{key_prefix}var"] = (
                    statistics.variance(values["var"]) if len(values["var"]) > 1 else 0
                )
            elif stats == "range":
                rhealpix_feature["properties"][f"{key_prefix}range"] = (
                    max(values["range"]) - min(values["range"])
                    if values["range"]
                    else 0
                )

            elif stats == "minority":
                freq = Counter(values["values"])
                min_item = min(freq.items(), key=lambda x: x[1])[0] if freq else None
                rhealpix_feature["properties"][f"{key_prefix}minority"] = min_item

            elif stats == "majority":
                freq = Counter(values["values"])
                max_item = max(freq.items(), key=lambda x: x[1])[0] if freq else None
                rhealpix_feature["properties"][f"{key_prefix}majority"] = max_item

            elif stats == "variety":
                rhealpix_feature["properties"][f"{key_prefix}variety"] = len(
                    set(values["values"])
                )

        rhealpix_features.append(rhealpix_feature)

    return rhealpix_features


def main():
    E = WGS84_ELLIPSOID
    rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)

    parser = argparse.ArgumentParser(description="Binning point to rHEALpix DGGS")
    parser.add_argument(
        "-point",
        "--point",
        type=str,
        required=True,
        help="GeoJSON file path (Point or MultiPoint)",
    )
    parser.add_argument(
        "-r", "--resolution", type=int, default=8, help="Resolution of the grid [0..15]"
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

    if resolution < 0 or resolution > 15:
        print("Error: Please select a resolution in [0..15].")
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
    rhealpix_features = rhealpix_bin(
        rhealpix_dggs, point_features, resolution, stats, category, field_name
    )

    out_name = os.path.splitext(os.path.basename(point))[0]
    out_path = f"{out_name}_bin_rhealpix_{resolution}_{stats}.geojson"

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(
            {"type": "FeatureCollection", "features": rhealpix_features}, f, indent=2
        )

    print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()

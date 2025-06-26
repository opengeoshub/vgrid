import platform
import argparse
import os
import json
import statistics
from collections import defaultdict, Counter
from shapely.geometry import Point, Polygon
from shapely.wkt import loads

from tqdm import tqdm
from vgrid.binning.bin_helper import get_default_stats_structure, append_stats_value
from vgrid.generator.settings import geodesic_dggs_to_feature

if platform.system() == "Windows":
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.generator.isea4tgrid import fix_isea4t_wkt, fix_isea4t_antimeridian_cells
    from vgrid.conversion.latlon2dggs import latlon2isea4t


def isea4t_bin(isea4t_dggs, point_features, resolution, stats, category, field_name):
    if platform.system() == "Windows":
        isea4t_bins = defaultdict(lambda: defaultdict(get_default_stats_structure))

        for point_feature in tqdm(point_features, desc="Binning points"):
            geom = point_feature["geometry"]
            props = point_feature.get("properties", {})

            if geom["type"] == "Point":
                point = Point(geom["coordinates"])
                isea4t_id = latlon2isea4t(point.y, point.x, resolution)
                append_stats_value(
                    isea4t_bins, isea4t_id, props, stats, category, field_name
                )

            elif geom["type"] == "MultiPoint":
                for coords in geom["coordinates"]:
                    point = Point(coords)
                    isea4t_id = latlon2isea4t(point.y, point.x, resolution)
                    append_stats_value(
                        isea4t_bins, isea4t_id, props, stats, category, field_name
                    )

        isea4t_features = []
        for isea4t_id, categories in isea4t_bins.items():
            cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(
                DggsCell(isea4t_id), ShapeStringFormat.WKT
            )
            cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
            if (
                isea4t_id.startswith("00")
                or isea4t_id.startswith("09")
                or isea4t_id.startswith("14")
                or isea4t_id.startswith("04")
                or isea4t_id.startswith("19")
            ):
                cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)

            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
            num_edges = 3

            if not cell_polygon.is_valid:
                continue

            isea4t_feature = geodesic_dggs_to_feature(
                "isea4t", isea4t_id, resolution, cell_polygon, num_edges
            )

            for cat, values in categories.items():
                key_prefix = "" if category is None else f"{cat}_"

                if stats == "count":
                    isea4t_feature["properties"][f"{key_prefix}count"] = values["count"]
                elif stats == "sum":
                    isea4t_feature["properties"][f"{key_prefix}sum"] = sum(
                        values["sum"]
                    )
                elif stats == "min":
                    isea4t_feature["properties"][f"{key_prefix}min"] = min(
                        values["min"]
                    )
                elif stats == "max":
                    isea4t_feature["properties"][f"{key_prefix}max"] = max(
                        values["max"]
                    )
                elif stats == "mean":
                    isea4t_feature["properties"][f"{key_prefix}mean"] = statistics.mean(
                        values["mean"]
                    )
                elif stats == "median":
                    isea4t_feature["properties"][f"{key_prefix}median"] = (
                        statistics.median(values["median"])
                    )
                elif stats == "std":
                    isea4t_feature["properties"][f"{key_prefix}std"] = (
                        statistics.stdev(values["std"]) if len(values["std"]) > 1 else 0
                    )
                elif stats == "var":
                    isea4t_feature["properties"][f"{key_prefix}var"] = (
                        statistics.variance(values["var"])
                        if len(values["var"]) > 1
                        else 0
                    )
                elif stats == "range":
                    isea4t_feature["properties"][f"{key_prefix}range"] = (
                        max(values["range"]) - min(values["range"])
                        if values["range"]
                        else 0
                    )

                elif stats == "minority":
                    freq = Counter(values["values"])
                    min_item = (
                        min(freq.items(), key=lambda x: x[1])[0] if freq else None
                    )
                    isea4t_feature["properties"][f"{key_prefix}minority"] = min_item

                elif stats == "majority":
                    freq = Counter(values["values"])
                    max_item = (
                        max(freq.items(), key=lambda x: x[1])[0] if freq else None
                    )
                    isea4t_feature["properties"][f"{key_prefix}majority"] = max_item

                elif stats == "variety":
                    isea4t_feature["properties"][f"{key_prefix}variety"] = len(
                        set(values["values"])
                    )

            isea4t_features.append(isea4t_feature)

        return isea4t_features


def main():
    parser = argparse.ArgumentParser(description="Binning point to isea4t DGGS")
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
        help="Resolution of the grid [0..25]",
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

    if platform.system() == "Windows":
        isea4t_dggs = Eaggr(Model.ISEA4T)
        args = parser.parse_args()

        resolution = args.resolution
        point = args.point
        stats = args.statistics
        category = args.category
        field_name = args.field

        if resolution < 0 or resolution > 25:
            print("Error: Please select a resolution in [0..25].")
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
        isea4t_features = isea4t_bin(
            isea4t_dggs, point_features, resolution, stats, category, field_name
        )

        out_name = os.path.splitext(os.path.basename(point))[0]
        out_path = f"{out_name}_bin_isea4t_{resolution}_{stats}.geojson"

        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(
                {"type": "FeatureCollection", "features": isea4t_features}, f, indent=2
            )

        print(f"GeoJSON saved as {out_path}")


if __name__ == "__main__":
    main()

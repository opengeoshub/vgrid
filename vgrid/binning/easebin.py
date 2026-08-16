"""
EASE Grid Binning Module

Bins point data into EASE (Equal-Area Scalable Earth) grid cells and computes
various statistics using an equal-area projection grid system.

Key Functions:
- ease_bin(): Core binning function with spatial joins and aggregation
- easebin(): Main user-facing function with multiple input/output formats
- easebin_cli(): Command-line interface for binning functionality
"""

import argparse
import geopandas as gpd
from vgrid.conversion.latlon2dggs import latlon2ease
from vgrid.conversion.dggs2geo.ease2geo import ease2geo
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_ease_resolution,
    aggregate_joined,
)
from vgrid.utils.constants import STATS_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS


def ease_bin(
    data,
    resolution,
    stats="count",
    category_col=None,
    numeric_col=None,
    lat_col="lat",
    lon_col="lon",
    **kwargs,
):
    """
    Bin point data into EASE grid cells and aggregate with pandas groupby.
    Supports custom stats (range, variety, minority, majority). Only
    Point/MultiPoint geometries are considered.

    Points are assigned to cells via latlon2ease (same as point2ease), then
    aggregated and joined to cell geometries from ease2geo.
    """
    resolution = validate_ease_resolution(resolution)
    points_gdf = process_input_data_bin(
        data, lat_col=lat_col, lon_col=lon_col, **kwargs
    )
    # Keep only points and multipoints; ignore others
    if not points_gdf.empty:
        points_gdf = points_gdf[
            points_gdf.geometry.geom_type.isin(["Point", "MultiPoint"])
        ].copy()
        if "MultiPoint" in set(points_gdf.geometry.geom_type.unique()):
            points_gdf = points_gdf.explode(index_parts=False, ignore_index=True)

    id_col = "ease"
    if points_gdf.empty:
        return gpd.GeoDataFrame(columns=[id_col, "geometry"], crs="EPSG:4326")

    points_gdf = points_gdf.copy()
    points_gdf[id_col] = [
        latlon2ease(geom.y, geom.x, resolution) for geom in points_gdf.geometry
    ]

    join_cols = []
    if category_col and category_col in points_gdf.columns:
        join_cols.append(category_col)
    if stats != "count" and numeric_col:
        if numeric_col not in points_gdf.columns:
            raise ValueError(f"numeric_col '{numeric_col}' not found in input data")
        join_cols.append(numeric_col)
    joined = points_gdf[[c for c in [id_col, *join_cols] if c is not None]]

    grouped = aggregate_joined(
        joined, id_col, stats=stats, category_col=category_col, numeric_col=numeric_col
    )
    grouped = grouped.reset_index()

    ease_rows = []
    for ease_id in grouped[id_col]:
        cell_polygon = ease2geo(ease_id)
        row = geodesic_dggs_to_geoseries(
            "ease", ease_id, resolution, cell_polygon, num_edges=4
        )
        ease_rows.append(row)
    grid_gdf = gpd.GeoDataFrame(ease_rows, geometry="geometry", crs="EPSG:4326")

    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(
        out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326"
    )
    return result_gdf


def easebin(
    data,
    resolution,
    stats="count",
    category_col=None,
    numeric_col=None,
    output_format="gpd",
    **kwargs,
):
    resolution = validate_ease_resolution(resolution)
    if stats != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")
    result_gdf = ease_bin(data, resolution, stats, category_col, numeric_col, **kwargs)
    output_name = None
    if output_format in OUTPUT_FORMATS:
        import os

        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_easebin_{resolution}"
        else:
            output_name = f"easebin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def easebin_cli():
    parser = argparse.ArgumentParser(description="Binning point data to EASE DGGS")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input data: GeoJSON file path, URL, or other vector file formats",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=4,
        help="Resolution of the grid [0..6]",
    )
    parser.add_argument(
        "-stats",
        "--statistics",
        choices=STATS_OPTIONS,
        default="count",
        help="Statistic option",
    )
    parser.add_argument(
        "-category",
        "--category",
        dest="category_col",
        required=False,
        help="Optional category field for grouping",
    )
    parser.add_argument(
        "-numeric_col",
        "--numeric_col",
        dest="numeric_col",
        required=False,
        help="Numeric field to compute statistics (required if stats != 'count')",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        required=False,
        default="gpd",
        choices=OUTPUT_FORMATS,
    )
    args = parser.parse_args()
    try:
        result = easebin(
            data=args.input,
            resolution=args.resolution,
            stats=args.statistics,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}")
        return


def main():
    easebin_cli()


if __name__ == "__main__":
    main()

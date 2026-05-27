"""
DGGRID Grid Binning Module

Bins point data into DGGRID cells and computes various statistics for
DGGRID types (e.g., ISEA7H, ISEA4T, FULLER4D, IGEO7).

Key Functions:
- dggrid_bin(): Core binning function with spatial joins and aggregation
- dggridbin(): Main user-facing function with multiple input/output formats
- dggridbin_cli(): Command-line interface for binning functionality
"""

import os
import argparse
import geopandas as gpd
from vgrid.generator.dggridgen import dggridgen
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_dggrid_type,
    validate_dggrid_resolution,
    create_dggrid_instance,
    aggregate_joined,
)
from vgrid.utils.constants import (
    STATS_OPTIONS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    DGGRID_TYPES,
)


def dggrid_bin(
    dggrid_instance,
    dggs_type: str,
    data,
    resolution: int,
    stats: str = "count",
    category_col: str | None = None,
    numeric_col: str | None = None,
    lat_col: str = "lat",
    lon_col: str = "lon",
    split_antimeridian: bool = False,
    aggregate: bool = False,
    **kwargs,
):
    """
    Bin point data into DGGRID cells and compute statistics.

    Parameters
    ----------
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
    aggregate : bool, optional
        When True (with split_antimeridian), dissolve split cell parts by
        global_id. Passed to dggridgen / generate_grid.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))

    if stats != "count" and not numeric_col:
        raise ValueError(
            "A numeric_col is required for statistics other than 'count'"
        )

    points_gdf = process_input_data_bin(
        data, lat_col=lat_col, lon_col=lon_col, **kwargs
    )
    if not points_gdf.empty:
        points_gdf = points_gdf[
            points_gdf.geometry.geom_type.isin(["Point", "MultiPoint"])
        ].copy()
        if "MultiPoint" in set(points_gdf.geometry.geom_type.unique()):
            points_gdf = points_gdf.explode(index_parts=False, ignore_index=True)

    minx, miny, maxx, maxy = points_gdf.total_bounds
    bbox = (minx, miny, maxx, maxy)
    id_col = "global_id"
    grid_gdf = dggridgen(
        dggrid_instance=dggrid_instance,
        dggs_type=dggs_type,
        resolution=resolution,
        output_format="gpd",
        bbox=bbox,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
    )
    if grid_gdf.crs is None:
        grid_gdf = grid_gdf.set_crs(points_gdf.crs)
    elif points_gdf.crs is not None and grid_gdf.crs != points_gdf.crs:
        grid_gdf = grid_gdf.to_crs(points_gdf.crs)

    join_cols = []
    if category_col and category_col in points_gdf.columns:
        join_cols.append(category_col)
    if stats != "count" and numeric_col:
        if numeric_col not in points_gdf.columns:
            raise ValueError(f"numeric_col '{numeric_col}' not found in input data")
        join_cols.append(numeric_col)
    left = points_gdf[[c for c in ["geometry", *join_cols] if c is not None]]
    joined = gpd.sjoin(
        left, grid_gdf[[id_col, "geometry"]], how="inner", predicate="within"
    )

    grouped = aggregate_joined(
        joined, id_col, stats=stats, category_col=category_col, numeric_col=numeric_col
    )
    grouped = grouped.reset_index()

    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    out = out.rename(columns={id_col: f"dggrid_{dggs_type.lower()}"})
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326")
    return result_gdf


def dggridbin(
    dggrid_instance,
    dggs_type: str,
    data,
    resolution: int,
    stats: str = "count",
    category_col: str | None = None,
    numeric_col: str | None = None,
    output_format: str = "gpd",
    split_antimeridian: bool = False,
    aggregate: bool = False,
    **kwargs,
):
    """
    Bin point data into DGGRID cells and compute statistics from various input formats.
    """
    result_gdf = dggrid_bin(
        dggrid_instance=dggrid_instance,
        dggs_type=dggs_type,
        data=data,
        resolution=resolution,
        stats=stats,
        category_col=category_col,
        numeric_col=numeric_col,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        **kwargs,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_{dggs_type.lower()}bin_{resolution}"
        else:
            output_name = f"{dggs_type.lower()}bin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def dggridbin_cli():
    """Command-line interface for DGGRID binning."""
    parser = argparse.ArgumentParser(description="Binning point data to DGGRID DGGS")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input data: GeoJSON file path, URL, or other vector file formats",
    )
    parser.add_argument(
        "-t",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGRID_TYPES.keys(),
        help="DGGRID type",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Resolution (integer)",
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
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )
    parser.add_argument(
        "-aggregate",
        "--aggregate",
        action="store_true",
        help="Aggregate the resulting polygons (dissolve by global_id when split_antimeridian is set)",
    )

    args = parser.parse_args()

    try:
        result = dggridbin(
            dggrid_instance=create_dggrid_instance(),
            dggs_type=args.dggs_type,
            data=args.input,
            resolution=args.resolution,
            stats=args.statistics,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
            split_antimeridian=args.split_antimeridian,
            aggregate=args.aggregate,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    dggridbin_cli()

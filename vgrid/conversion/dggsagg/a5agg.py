"""
A5 Aggregate Module

Roll A5 cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    a5_agg: Group A5 hex IDs by parent cell at a resolution
    a5agg: Read cells into a GeoDataFrame and aggregate by parent
    a5agg_cli: Command-line interface
"""

import os
import json
import argparse
from collections import defaultdict

import geopandas as gpd
import pandas as pd
import a5
from tqdm import tqdm
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    convert_to_output_format,
    prepare_compact_bags,
    validate_a5_resolution,
)
from vgrid.utils.constants import (
    AGG_OPTIONS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)
from vgrid.conversion.dggs2geo.a52geo import a52geo


def a5_parent_at(a5_hex, resolution):
    """Parent of ``a5_hex`` at ``resolution``, or the cell itself when already there."""
    u = a5.hex_to_u64(str(a5_hex))
    cell_res = a5.get_resolution(u)
    if cell_res < resolution:
        raise ValueError(
            f"A5 cell {a5_hex} is resolution {cell_res}, coarser than "
            f"parent resolution {resolution}."
        )
    if cell_res == resolution:
        return a5.u64_to_hex(u)
    return a5.u64_to_hex(a5.cell_to_parent(u, resolution))


def a5_agg(a5_hexes, resolution, bags=None, verbose=True):
    """Group A5 cell hex IDs by their parent at ``resolution``.

    Every cell is assigned to ``a5.cell_to_parent(cell, resolution)``, whether
    or not its siblings are present. Cells already at ``resolution`` stay put.
    Cells coarser than ``resolution`` raise ``ValueError``.

    Parameters
    ----------
    a5_hexes : list of str
        A5 cell hex IDs. Mixed resolutions finer than ``resolution`` are allowed.
    resolution : int
        Parent A5 resolution to group by.
    bags : dict of list, optional
        Per-cell lists of original values. Child lists are concatenated onto
        the parent. Mutated so remaining keys are the parent IDs.
    verbose : bool, default True
        Show a tqdm progress bar.

    Returns
    -------
    list of str
        Sorted parent A5 cell hex IDs.
    """
    resolution = validate_a5_resolution(resolution)
    parent_bags = defaultdict(list) if bags is not None else None
    parents = set()
    for a5_hex in tqdm(
        a5_hexes,
        desc="Aggregating A5",
        unit=" cells",
        disable=not verbose,
    ):
        parent = a5_parent_at(a5_hex, resolution)
        parents.add(parent)
        if parent_bags is not None:
            parent_bags[parent].extend(bags.get(a5_hex, []))
    if bags is not None:
        bags.clear()
        bags.update(parent_bags)
    return sorted(parents)


def a5agg(
    input_data,
    resolution,
    a5_hex="a5",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    options=None,
    split_antimeridian=False,
    cell_metrics=False,
    verbose=True,
):
    """Aggregate A5 cells by parent at ``resolution``.

    Reads the input into a GeoDataFrame, maps each cell to
    ``a5.cell_to_parent(cell, resolution)``, and combines values with ``agg``.
    Missing siblings are still rolled up.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing A5 cell hex IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of A5 cell hex IDs
    resolution : int
        Parent A5 resolution to group by. Must be between 0 and the finest
        input cell.
    a5_hex : str, optional
        Name of the column containing A5 hex IDs. Defaults to "a5".
    agg : str, default "count"
        Aggregation applied to original values in each parent. Same options as
        ``a5compact`` (``count``, ``min``, ``max``, ``sum``, ``mean``,
        ``median``, ``std``, ``var``, ``range``, ``minority``, ``majority``,
        ``variety``).
    numeric_col : str, optional
        Numeric field to aggregate. Required when ``agg`` is not ``"count"``;
        ignored when ``agg`` is ``"count"``.
    options : dict, optional
        Options for a52geo. Used when ``cell_metrics`` is true.
    split_antimeridian : bool, default False
        Split parent polygons at the antimeridian. Used when ``cell_metrics``
        is true.
    cell_metrics : bool, default False
        When true, each parent row includes cell geometry and geodesic metrics
        from ``geodesic_dggs_to_geoseries``. Otherwise the result has only the
        A5 id column and the aggregated column.
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        Parent A5 cells in the specified format, or None if no valid cells found.
    """
    if a5_hex is None:
        a5_hex = "a5"
    resolution = validate_a5_resolution(resolution)
    bags, agg_col = prepare_compact_bags(
        input_data,
        a5_hex,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="A5 cells",
    )
    if bags is None:
        print(f"No A5 IDs found in <{a5_hex}> field.")
        return

    parent_hexes = a5_agg(list(bags.keys()), resolution, bags=bags, verbose=verbose)
    if not parent_hexes:
        return None

    num_edges = 3 if resolution == 1 else 5
    rows = []
    for parent_hex in tqdm(
        parent_hexes,
        desc="Building A5 aggregate",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            if cell_metrics:
                cell_polygon = a52geo(
                    parent_hex, options, split_antimeridian=split_antimeridian
                )
                row = geodesic_dggs_to_geoseries(
                    "a5", parent_hex, resolution, cell_polygon, num_edges
                )
            else:
                row = {a5_hex: parent_hex}
            row[agg_col] = aggregate_values(bags.get(parent_hex, []), agg)
            rows.append(row)
        except Exception:
            continue
    if cell_metrics:
        out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    else:
        out_gdf = pd.DataFrame(rows)

    if not cell_metrics and output_format in (
        "gpd",
        "geopandas",
        "gdf",
        "geodataframe",
    ):
        return out_gdf

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_a5_agg"
        else:
            output_name = "a5_agg"

    return convert_to_output_format(out_gdf, output_format, output_name)


def a5agg_cli():
    """Command-line interface for a5agg."""
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate A5 cells by parent resolution. "
            "Unlike a5compact, incomplete child sets are still rolled up."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input A5 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Parent A5 resolution to group by (a5.cell_to_parent).",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="A5 Hex field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable Antimeridian splitting",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
        "Example: '{\"segments\": 1000}'",
    )
    parser.add_argument(
        "-agg",
        "--agg",
        choices=AGG_OPTIONS,
        default="count",
        help="Aggregation option",
    )
    parser.add_argument(
        "-numeric_col",
        "--numeric_col",
        dest="numeric_col",
        required=False,
        help="Numeric field to aggregate (required if agg != 'count')",
    )
    parser.add_argument(
        "-cell_metrics",
        "--cell_metrics",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Include cell geometry and geodesic metrics. Default: omit them.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress bar (default: True). Use --no-verbose to hide it.",
    )

    args = parser.parse_args()
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error parsing options JSON: {e}")
            return

    result = a5agg(
        args.input,
        resolution=args.resolution,
        a5_hex=args.cellid,
        output_format=args.output_format,
        options=options,
        split_antimeridian=args.split_antimeridian,
        agg=args.agg,
        numeric_col=args.numeric_col,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)

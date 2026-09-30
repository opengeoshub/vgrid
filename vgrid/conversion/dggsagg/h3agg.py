"""
H3 Aggregate Module

Roll H3 cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    h3_agg: Group H3 IDs by parent cell at a resolution
    h3agg: Read cells into a GeoDataFrame and aggregate by parent
    h3agg_cli: Command-line interface
"""

import os
import argparse
from collections import defaultdict

import geopandas as gpd
import pandas as pd
import h3
from tqdm import tqdm
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    convert_to_output_format,
    prepare_compact_bags,
    validate_h3_resolution,
)
from vgrid.utils.constants import (
    AGG_OPTIONS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    FIX_ANTIMERIDIAN_CHOICES,
)
from vgrid.conversion.dggs2geo.h32geo import h32geo


def h3_parent_at(h3_id, resolution):
    """Parent of ``h3_id`` at ``resolution``, or the cell itself when already there."""
    cell_res = h3.get_resolution(h3_id)
    if cell_res < resolution:
        raise ValueError(
            f"H3 cell {h3_id} is resolution {cell_res}, coarser than "
            f"parent resolution {resolution}."
        )
    if cell_res == resolution:
        return h3_id
    return h3.cell_to_parent(h3_id, resolution)


def h3_agg(h3_ids, resolution, bags=None, verbose=True):
    """Group H3 cell IDs by their parent at ``resolution``.

    Every cell is assigned to ``h3.cell_to_parent(cell, resolution)``, whether
    or not its siblings are present. Cells already at ``resolution`` stay put.
    Cells coarser than ``resolution`` raise ``ValueError``.

    Parameters
    ----------
    h3_ids : list of str
        H3 cell IDs. Mixed resolutions finer than ``resolution`` are allowed.
    resolution : int
        Parent H3 resolution to group by.
    bags : dict of list, optional
        Per-cell lists of original values. Child lists are concatenated onto
        the parent. Mutated so remaining keys are the parent IDs.
    verbose : bool, default True
        Show a tqdm progress bar.

    Returns
    -------
    list of str
        Sorted parent H3 cell IDs.
    """
    resolution = validate_h3_resolution(resolution)
    parent_bags = defaultdict(list) if bags is not None else None
    parents = set()
    for h3_id in tqdm(
        h3_ids,
        desc="Aggregating H3",
        unit=" cells",
        disable=not verbose,
    ):
        parent = h3_parent_at(h3_id, resolution)
        parents.add(parent)
        if parent_bags is not None:
            parent_bags[parent].extend(bags.get(h3_id, []))
    if bags is not None:
        bags.clear()
        bags.update(parent_bags)
    return sorted(parents)


def h3agg(
    input_data,
    resolution,
    h3_id="h3",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian="shift",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate H3 cells by parent at ``resolution``.

    Reads the input into a GeoDataFrame, maps each cell to
    ``h3.cell_to_parent(cell, resolution)``, and combines values with ``agg``.
    Missing siblings are still rolled up. This matches::

        SELECT h3_cell_to_parent(cell, resolution) AS cell,
               SUM(value) AS value
        FROM cells
        GROUP BY 1

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing H3 cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of H3 cell IDs
    resolution : int
        Parent H3 resolution to group by. Must be between 0 and the finest
        input cell.
    h3_id : str, optional
        Name of the column containing H3 cell IDs. Defaults to "h3".
    agg : str, default "count"
        Aggregation applied to original values in each parent. Same options as
        ``h3compact`` (``count``, ``min``, ``max``, ``sum``, ``mean``,
        ``median``, ``std``, ``var``, ``range``, ``minority``, ``majority``,
        ``variety``).
    numeric_col : str, optional
        Numeric field to aggregate. Required when ``agg`` is not ``"count"``;
        ignored when ``agg`` is ``"count"``.
    fix_antimeridian : str, default "shift"
        Antimeridian fixing for parent cells: shift, shift_balanced, shift_west,
        shift_east, split, none. Used when ``cell_metrics`` is true.
    cell_metrics : bool, default False
        When true, each parent row includes cell geometry and geodesic metrics
        from ``geodesic_dggs_to_geoseries``. Otherwise the result has only the
        H3 id column and the aggregated column.
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
        Parent H3 cells in the specified format, or None if no valid cells found.
    """
    if h3_id is None:
        h3_id = "h3"
    resolution = validate_h3_resolution(resolution)
    bags, agg_col = prepare_compact_bags(
        input_data,
        h3_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="H3 cells",
    )
    if bags is None:
        print(f"No H3 IDs found in <{h3_id}> field.")
        return

    parent_ids = h3_agg(list(bags.keys()), resolution, bags=bags, verbose=verbose)
    if not parent_ids:
        return None

    rows = []
    for parent_id in tqdm(
        parent_ids,
        desc="Building H3 aggregate",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            if cell_metrics:
                cell_polygon = h32geo(parent_id, fix_antimeridian=fix_antimeridian)
                num_edges = 5 if h3.is_pentagon(parent_id) else 6
                row = geodesic_dggs_to_geoseries(
                    "h3", parent_id, resolution, cell_polygon, num_edges
                )
            else:
                row = {h3_id: parent_id}
            row[agg_col] = aggregate_values(bags.get(parent_id, []), agg)
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
            output_name = f"{base}_h3_agg"
        else:
            output_name = "h3_agg"

    return convert_to_output_format(out_gdf, output_format, output_name)


def h3agg_cli():
    """Command-line interface for h3agg."""
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate H3 cells by parent resolution. "
            "Unlike h3compact, incomplete child sets are still rolled up."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input H3 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Parent H3 resolution to group by (h3.cell_to_parent).",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="H3 ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default="shift",
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none. Default: shift.",
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
    result = h3agg(
        args.input,
        resolution=args.resolution,
        h3_id=args.cellid,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        agg=args.agg,
        numeric_col=args.numeric_col,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)

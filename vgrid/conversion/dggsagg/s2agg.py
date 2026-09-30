"""
S2 Aggregate Module

Roll S2 cells up to a parent resolution and aggregate values there.
Incomplete child sets are included. This is a group-by parent, not compaction.

Key Functions:
    s2_agg: Group S2 tokens by parent cell at a resolution
    s2agg: Read cells into a GeoDataFrame and aggregate by parent
    s2agg_cli: Command-line interface
"""

import os
import argparse
from collections import defaultdict

import geopandas as gpd
import pandas as pd
from tqdm import tqdm
from vgrid.dggs import s2
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    convert_to_output_format,
    prepare_compact_bags,
    validate_s2_resolution,
)
from vgrid.utils.constants import (
    AGG_OPTIONS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    FIX_ANTIMERIDIAN_CHOICES,
)
from vgrid.conversion.dggs2geo.s22geo import s22geo


def s2_parent_at(s2_token, resolution):
    """Parent of ``s2_token`` at ``resolution``, or the cell itself when already there."""
    cell_id = s2.CellId.from_token(str(s2_token))
    cell_res = cell_id.level()
    if cell_res < resolution:
        raise ValueError(
            f"S2 cell {s2_token} is resolution {cell_res}, coarser than "
            f"parent resolution {resolution}."
        )
    if cell_res == resolution:
        return cell_id.to_token()
    return cell_id.parent(resolution).to_token()


def s2_agg(s2_tokens, resolution, bags=None, verbose=True):
    """Group S2 cell tokens by their parent at ``resolution``.

    Every cell is assigned to ``CellId.parent(resolution)``, whether or not its
    siblings are present. Cells already at ``resolution`` stay put. Cells
    coarser than ``resolution`` raise ``ValueError``.

    Parameters
    ----------
    s2_tokens : list of str
        S2 cell tokens. Mixed resolutions finer than ``resolution`` are allowed.
    resolution : int
        Parent S2 resolution (level) to group by.
    bags : dict of list, optional
        Per-cell lists of original values. Child lists are concatenated onto
        the parent. Mutated so remaining keys are the parent tokens.
    verbose : bool, default True
        Show a tqdm progress bar.

    Returns
    -------
    list of str
        Sorted parent S2 cell tokens.
    """
    resolution = validate_s2_resolution(resolution)
    parent_bags = defaultdict(list) if bags is not None else None
    parents = set()
    for s2_token in tqdm(
        s2_tokens,
        desc="Aggregating S2",
        unit=" cells",
        disable=not verbose,
    ):
        parent = s2_parent_at(s2_token, resolution)
        parents.add(parent)
        if parent_bags is not None:
            parent_bags[parent].extend(bags.get(s2_token, []))
    if bags is not None:
        bags.clear()
        bags.update(parent_bags)
    return sorted(parents)


def s2agg(
    input_data,
    resolution,
    s2_token="s2",
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian="shift",
    cell_metrics=False,
    verbose=True,
):
    """Aggregate S2 cells by parent at ``resolution``.

    Reads the input into a GeoDataFrame, maps each cell to
    ``CellId.parent(resolution)``, and combines values with ``agg``.
    Missing siblings are still rolled up.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing S2 cell tokens. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of S2 cell tokens
    resolution : int
        Parent S2 resolution to group by. Must be between 0 and the finest
        input cell.
    s2_token : str, optional
        Name of the column containing S2 cell tokens. Defaults to "s2".
    agg : str, default "count"
        Aggregation applied to original values in each parent. Same options as
        ``s2compact`` (``count``, ``min``, ``max``, ``sum``, ``mean``,
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
        S2 token column and the aggregated column.
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
        Parent S2 cells in the specified format, or None if no valid cells found.
    """
    if s2_token is None:
        s2_token = "s2"
    resolution = validate_s2_resolution(resolution)
    bags, agg_col = prepare_compact_bags(
        input_data,
        s2_token,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="S2 cells",
    )
    if bags is None:
        print(f"No S2 tokens found in <{s2_token}> field.")
        return

    parent_tokens = s2_agg(list(bags.keys()), resolution, bags=bags, verbose=verbose)
    if not parent_tokens:
        return None

    rows = []
    for parent_token in tqdm(
        parent_tokens,
        desc="Building S2 aggregate",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            if cell_metrics:
                cell_polygon = s22geo(parent_token, fix_antimeridian=fix_antimeridian)
                row = geodesic_dggs_to_geoseries(
                    "s2", parent_token, resolution, cell_polygon, 4
                )
            else:
                row = {s2_token: parent_token}
            row[agg_col] = aggregate_values(bags.get(parent_token, []), agg)
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
            output_name = f"{base}_s2_agg"
        else:
            output_name = "s2_agg"

    return convert_to_output_format(out_gdf, output_format, output_name)


def s2agg_cli():
    """Command-line interface for s2agg."""
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate S2 cells by parent resolution. "
            "Unlike s2compact, incomplete child sets are still rolled up."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input S2 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Parent S2 resolution to group by (CellId.parent).",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="S2 token field")
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
    result = s2agg(
        args.input,
        resolution=args.resolution,
        s2_token=args.cellid,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        agg=args.agg,
        numeric_col=args.numeric_col,
        cell_metrics=args.cell_metrics,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)

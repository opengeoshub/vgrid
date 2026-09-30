"""
Shared helpers for the DGGS aggregate modules.

Each ``<dggs>agg`` module supplies how to find a cell's parent at a target
resolution and how to build a parent row with geometry and metrics; the
grouping, aggregation, output, and CLI plumbing live here.
"""

import os
import argparse
from collections import defaultdict

import geopandas as gpd
import pandas as pd
from tqdm import tqdm
from vgrid.utils.io import (
    aggregate_values,
    convert_to_output_format,
    prepare_compact_bags,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS

_GDF_FORMATS = ("gpd", "geopandas", "gdf", "geodataframe")


def climb_to_parent(cell_id, resolution, get_resolution, parent_fn, label):
    """Walk ``parent_fn`` up from ``cell_id`` until ``resolution`` is reached.

    Raises ``ValueError`` when the cell is coarser than ``resolution`` or runs
    out of parents before reaching it.
    """
    cell_res = get_resolution(cell_id)
    if cell_res < resolution:
        raise ValueError(
            f"{label} cell {cell_id} is resolution {cell_res}, coarser than "
            f"parent resolution {resolution}."
        )
    current = cell_id
    while cell_res > resolution:
        current = parent_fn(current)
        if current is None:
            raise ValueError(
                f"{label} cell {cell_id} has no parent at resolution {resolution}."
            )
        cell_res = get_resolution(current)
    if cell_res != resolution:
        raise ValueError(
            f"{label} cell {cell_id} skips parent resolution {resolution}."
        )
    return current


def check_not_coarser(cell_id, cell_res, resolution, label):
    if cell_res < resolution:
        raise ValueError(
            f"{label} cell {cell_id} is resolution {cell_res}, coarser than "
            f"parent resolution {resolution}."
        )


def group_by_parent(cell_ids, parent_of, label, bags=None, verbose=True):
    """Group cells by ``parent_of(cell_id)``.

    ``bags`` holds per-cell lists of original values; child lists are
    concatenated onto the parent and ``bags`` is mutated so its keys are the
    parent IDs. Returns the sorted parent IDs.
    """
    parent_bags = defaultdict(list) if bags is not None else None
    parents = set()
    for cell_id in tqdm(
        cell_ids,
        desc=f"Aggregating {label}",
        unit=" cells",
        disable=not verbose,
    ):
        parent = parent_of(cell_id)
        parents.add(parent)
        if parent_bags is not None:
            parent_bags[parent].extend(bags.get(cell_id, []))
    if bags is not None:
        bags.clear()
        bags.update(parent_bags)
    return sorted(parents)


def build_agg_rows(
    parent_ids,
    bags,
    agg,
    agg_col,
    id_col,
    row_fn,
    cell_metrics,
    label,
    verbose=True,
):
    """One row per parent: ``row_fn(parent)`` when ``cell_metrics`` else the id."""
    rows = []
    for parent_id in tqdm(
        parent_ids,
        desc=f"Building {label} aggregate",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            row = row_fn(parent_id) if cell_metrics else {id_col: parent_id}
            row[agg_col] = aggregate_values(bags.get(parent_id, []), agg)
            rows.append(row)
        except Exception:
            continue
    return rows


def finish_agg_output(rows, cell_metrics, output_format, input_data, name):
    """GeoDataFrame with metrics, plain DataFrame without, then convert."""
    if cell_metrics:
        out_df = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    else:
        out_df = pd.DataFrame(rows)

    if not cell_metrics and output_format in _GDF_FORMATS:
        return out_df

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_{name}_agg"
        else:
            output_name = f"{name}_agg"
    return convert_to_output_format(out_df, output_format, output_name)


def run_agg(
    input_data,
    id_col,
    label,
    name,
    parent_of,
    row_fn,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    cell_metrics=False,
    verbose=True,
):
    """Read cells, group them by ``parent_of``, aggregate, and convert."""
    bags, agg_col = prepare_compact_bags(
        input_data,
        id_col,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label=f"{label} cells",
    )
    if bags is None:
        print(f"No {label} IDs found in <{id_col}> field.")
        return None

    parent_ids = group_by_parent(
        list(bags.keys()), parent_of, label, bags=bags, verbose=verbose
    )
    if not parent_ids:
        return None

    rows = build_agg_rows(
        parent_ids,
        bags,
        agg,
        agg_col,
        id_col,
        row_fn,
        cell_metrics,
        label,
        verbose=verbose,
    )
    return finish_agg_output(rows, cell_metrics, output_format, input_data, name)


def agg_arg_parser(label, compact_cmd):
    """Parser with the options every ``<dggs>agg`` CLI shares."""
    parser = argparse.ArgumentParser(
        description=(
            f"Aggregate {label} cells by parent resolution. "
            f"Unlike {compact_cmd}, incomplete child sets are still rolled up."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help=f"Input {label} (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help=f"Parent {label} resolution to group by.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help=f"{label} ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
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
        help="Include cell geometry and metrics. Default: omit them.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress bar (default: True). Use --no-verbose to hide it.",
    )
    return parser


def print_structured(output_format, result):
    if output_format in STRUCTURED_FORMATS:
        print(result)

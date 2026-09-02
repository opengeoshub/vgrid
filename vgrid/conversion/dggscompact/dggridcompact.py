"""
DGGRID Compact Module

This module provides functionality to compact and expand DGGRID cells
with flexible input and output formats.

Key Functions:
    dggridcompact: Compact a set of DGGRID cells to a coarser minimal covering set
    dggridexpand: Expand (uncompact) a set of DGGRID cells to a target resolution
    dggridcompact_cli: Command-line interface for compaction
    dggridexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import json
import geopandas as gpd
from tqdm import tqdm
from vgrid.generator.dggridgen import dggridgen
from vgrid.conversion.dggs2geo.dggrid2geo import dggrid2geo
from vgrid.utils.geometry import dggrid_num_edges, geodesic_dggs_metrics
from vgrid.utils.io import (
    aggregate_values,
    convert_to_output_format,
    create_dggrid_instance,
    prepare_compact_bags,
    process_input_data_compact,
    validate_dggrid_resolution,
    validate_dggrid_type,
)
from vgrid.utils.constants import AGG_OPTIONS, DGGRID_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS


def _extract_geom(cell_gdf):
    if isinstance(cell_gdf, gpd.GeoDataFrame) and not cell_gdf.empty:
        return cell_gdf.iloc[0].geometry
    return None


def _resolve_cell_geometry(
    dggrid_instance,
    dggs_type,
    cell_id,
    preferred_resolution,
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    min_res = int(DGGRID_TYPES[dggs_type]["min_res"])
    max_res = int(DGGRID_TYPES[dggs_type]["max_res"])
    candidates = [preferred_resolution] + [
        r for r in range(min_res, max_res + 1) if r != preferred_resolution
    ]
    for res in candidates:
        try:
            cell_gdf = dggrid2geo(
                dggrid_instance,
                dggs_type,
                cell_id,
                res,
                split_antimeridian=split_antimeridian,
                aggregate=aggregate,
                options=options,
            )
            geom = _extract_geom(cell_gdf)
            if geom is not None:
                return geom, res
        except Exception:
            continue
    return None, None


def _cells_to_gdf(
    dggrid_instance,
    dggs_type,
    cell_ids,
    resolution,
    split_antimeridian=False,
    aggregate=False,
    options=None,
    verbose=True,
    desc="Building DGGRID cells",
):
    rows = []
    id_col = f"dggrid_{dggs_type.lower()}"
    for cell_id in tqdm(
        cell_ids, desc=desc, unit=" cells", disable=not verbose
    ):
        try:
            geom, cell_resolution = _resolve_cell_geometry(
                dggrid_instance,
                dggs_type,
                cell_id,
                preferred_resolution=resolution,
                split_antimeridian=split_antimeridian,
                aggregate=aggregate,
                options=options,
            )
            if geom is None or cell_resolution is None:
                continue
            centroid_lat, centroid_lon, avg_edge_len, cell_area, cell_perimeter = (
                geodesic_dggs_metrics(geom, dggrid_num_edges(dggs_type))
            )
            rows.append(
                {
                    id_col: cell_id,
                    "resolution": cell_resolution,
                    "center_lat": centroid_lat,
                    "center_lon": centroid_lon,
                    "avg_edge_len": avg_edge_len,
                    "cell_area": cell_area,
                    "cell_perimeter": cell_perimeter,
                    "geometry": geom,
                }
            )
        except Exception:
            continue
    return gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")


def dggrid_compact(
    dggrid_instance,
    dggs_type,
    cell_ids,
    resolution,
    depth=-1,
    bags=None,
    verbose=True,
):
    """
    Compact DGGRID cells by iteratively replacing complete child sets with parents.

    Parent/child mapping is spatial (generate grids + sjoin), so this path does
    not use ``compact_cells``. ``depth`` limits how many parent-level iterations
    run: ``0`` is a no-op, ``-1`` compact until min resolution or no change.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))
    min_res = int(DGGRID_TYPES[dggs_type]["min_res"])

    active_ids = set(str(v) for v in cell_ids)
    if depth == 0 or not active_ids:
        return sorted(active_ids)

    current_res = resolution
    max_iterations = None if depth < 0 else depth
    iterations = 0

    while current_res > min_res and active_ids:
        if max_iterations is not None and iterations >= max_iterations:
            break
        sample_gdf = _cells_to_gdf(
            dggrid_instance=dggrid_instance,
            dggs_type=dggs_type,
            cell_ids=list(active_ids),
            resolution=current_res,
            verbose=verbose,
            desc="Resolving DGGRID cells",
        )
        if sample_gdf.empty:
            break
        minx, miny, maxx, maxy = sample_gdf.total_bounds
        bbox = (minx, miny, maxx, maxy)

        parent_gdf = dggridgen(
            dggrid_instance,
            dggs_type,
            current_res - 1,
            bbox=bbox,
            output_format="gpd",
        )
        child_gdf = dggridgen(
            dggrid_instance,
            dggs_type,
            current_res,
            bbox=bbox,
            output_format="gpd",
        )
        if parent_gdf.empty or child_gdf.empty:
            break

        # Map candidate children to parents by centroid containment.
        child_cent = child_gdf.copy()
        child_cent["geometry"] = child_cent.geometry.representative_point()
        child_cent = child_cent.rename(columns={"global_id": "child_id"})
        parent_ids = parent_gdf[["global_id", "geometry"]].rename(
            columns={"global_id": "parent_id"}
        )
        mapping = gpd.sjoin(
            child_cent[["child_id", "geometry"]],
            parent_ids,
            how="inner",
            predicate="within",
        )
        if mapping.empty:
            break

        parent_to_all = mapping.groupby("parent_id")["child_id"].apply(set)
        child_universe = set(int(v) for v in child_gdf["global_id"].tolist())
        active_num = set()
        for v in active_ids:
            try:
                iv = int(v)
                if iv in child_universe:
                    active_num.add(iv)
            except Exception:
                continue
        if not active_num:
            break

        changed = False
        new_active = set(active_ids)
        for parent_id, all_children in parent_to_all.items():
            if all_children and all_children.issubset(active_num):
                # Replace all children by parent id
                new_active.difference_update(str(c) for c in all_children)
                new_active.add(str(parent_id))
                if bags is not None:
                    merged = []
                    parent_key = str(parent_id)
                    if parent_key in bags and parent_id not in all_children:
                        merged.extend(bags[parent_key])
                    for c in all_children:
                        merged.extend(bags.pop(str(c), []))
                    bags[parent_key] = merged
                changed = True
        if not changed:
            break
        active_ids = new_active
        current_res -= 1
        iterations += 1

    if bags is not None:
        for gone in list(bags.keys()):
            if gone not in active_ids:
                bags.pop(gone, None)

    return sorted(active_ids)


def dggrid_expand(
    dggrid_instance,
    dggs_type,
    cell_ids,
    resolution,
    target_resolution,
):
    """
    Expand DGGRID cells to target resolution using spatial containment.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))
    target_resolution = validate_dggrid_resolution(dggs_type, int(target_resolution))
    if target_resolution < resolution:
        raise ValueError(
            f"Target resolution ({target_resolution}) must be >= input resolution ({resolution})."
        )
    if target_resolution == resolution:
        return sorted(set(str(v) for v in cell_ids))

    # Build source geometry robustly, allowing mixed-resolution compacted IDs.
    src_gdf = _cells_to_gdf(
        dggrid_instance=dggrid_instance,
        dggs_type=dggs_type,
        cell_ids=cell_ids,
        resolution=resolution,
    )
    if src_gdf.empty:
        return []
    minx, miny, maxx, maxy = src_gdf.total_bounds
    bbox = (minx, miny, maxx, maxy)

    tgt_gdf = dggridgen(
        dggrid_instance,
        dggs_type,
        target_resolution,
        bbox=bbox,
        output_format="gpd",
    )
    if tgt_gdf.empty:
        return []

    if tgt_gdf.crs is None and src_gdf.crs is not None:
        tgt_gdf = tgt_gdf.set_crs(src_gdf.crs)
    elif src_gdf.crs is None and tgt_gdf.crs is not None:
        src_gdf = src_gdf.set_crs(tgt_gdf.crs)
    elif (
        tgt_gdf.crs is not None
        and src_gdf.crs is not None
        and tgt_gdf.crs != src_gdf.crs
    ):
        tgt_gdf = tgt_gdf.to_crs(src_gdf.crs)

    tgt_cent = tgt_gdf.copy()
    tgt_cent["geometry"] = tgt_cent.geometry.representative_point()
    joined = gpd.sjoin(
        tgt_cent[["global_id", "geometry"]],
        src_gdf[["geometry"]],
        how="inner",
        predicate="within",
    )
    if joined.empty:
        return []
    expanded = sorted(str(v) for v in joined["global_id"].drop_duplicates().tolist())
    return expanded


def dggridcompact(
    dggrid_instance,
    dggs_type,
    input_data,
    resolution,
    dggrid_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    split_antimeridian=False,
    aggregate=False,
    options=None,
    verbose=True,
):
    """
    Compact DGGRID cells to their covering set at a given parent depth.

    Parent/child mapping is spatial. When a complete child set is replaced by
    its parent, original child values are combined with ``agg``.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))
    if dggrid_id is None:
        dggrid_id = f"dggrid_{dggs_type.lower()}"

    bags, agg_col = prepare_compact_bags(
        input_data,
        dggrid_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="DGGRID cells",
    )
    if bags is None:
        print(f"No DGGRID IDs found in <{dggrid_id}> field.")
        return None

    compact_ids = dggrid_compact(
        dggrid_instance=dggrid_instance,
        dggs_type=dggs_type,
        cell_ids=list(bags.keys()),
        resolution=resolution,
        depth=depth,
        bags=bags,
        verbose=verbose,
    )
    if not compact_ids:
        return None

    out_gdf = _cells_to_gdf(
        dggrid_instance,
        dggs_type,
        compact_ids,
        resolution=resolution,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
        verbose=verbose,
        desc="Building DGGRID compact",
    )
    id_col = f"dggrid_{dggs_type.lower()}"
    if not out_gdf.empty and id_col in out_gdf.columns:
        out_gdf[agg_col] = [
            aggregate_values(bags.get(str(cid), []), agg) for cid in out_gdf[id_col]
        ]

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_dggrid_compacted"
        else:
            output_name = "dggrid_compacted"
    return convert_to_output_format(out_gdf, output_format, output_name)


def dggridexpand(
    dggrid_instance,
    dggs_type,
    input_data,
    resolution,
    target_resolution,
    dggrid_id=None,
    output_format="gpd",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Expand (uncompact) DGGRID cells to a target resolution.
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, int(resolution))
    target_resolution = validate_dggrid_resolution(dggs_type, int(target_resolution))
    if dggrid_id is None:
        dggrid_id = f"dggrid_{dggs_type.lower()}"

    gdf = process_input_data_compact(input_data, dggrid_id)
    cell_ids = gdf[dggrid_id].drop_duplicates().tolist()
    if not cell_ids:
        print(f"No DGGRID IDs found in <{dggrid_id}> field.")
        return None

    expanded_ids = dggrid_expand(
        dggrid_instance=dggrid_instance,
        dggs_type=dggs_type,
        cell_ids=cell_ids,
        resolution=resolution,
        target_resolution=target_resolution,
    )
    if not expanded_ids:
        return None

    out_gdf = _cells_to_gdf(
        dggrid_instance,
        dggs_type,
        expanded_ids,
        resolution=target_resolution,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_dggrid_expanded"
        else:
            output_name = "dggrid_expanded"
    return convert_to_output_format(out_gdf, output_format, output_name)


def dggridcompact_cli():
    parser = argparse.ArgumentParser(description="DGGRID Compact")
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Input DGGRID file"
    )
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGRID_TYPES.keys(),
        help="DGGRID type",
    )
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution"
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="DGGRID ID field")
    parser.add_argument(
        "-f", "--output_format", type=str, default="gpd", choices=OUTPUT_FORMATS
    )
    parser.add_argument(
        "-split", "--split_antimeridian", action="store_true", default=False
    )
    parser.add_argument("-aggregate", "--aggregate", action="store_true")
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string options for dggrid2geo/dggridgen",
    )
    parser.add_argument(
        "-d",
        "--depth",
        type=int,
        default=-1,
        help="Compaction depth: 0 = no-op, -1 = compact fully (default), "
        "1 = direct parent, 2 = grandparent, ...",
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
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    dggrid_instance = create_dggrid_instance()
    result = dggridcompact(
        dggrid_instance=dggrid_instance,
        dggs_type=args.dggs_type,
        input_data=args.input,
        resolution=args.resolution,
        dggrid_id=args.cellid,
        output_format=args.output_format,
        split_antimeridian=args.split_antimeridian,
        aggregate=args.aggregate,
        options=options,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


def dggridexpand_cli():
    parser = argparse.ArgumentParser(description="DGGRID Expand (Uncompact)")
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="Input DGGRID file"
    )
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGRID_TYPES.keys(),
        help="DGGRID type",
    )
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Input resolution"
    )
    parser.add_argument(
        "-tr", "--target_resolution", type=int, required=True, help="Target resolution"
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="DGGRID ID field")
    parser.add_argument(
        "-f", "--output_format", type=str, default="gpd", choices=OUTPUT_FORMATS
    )
    parser.add_argument(
        "-split", "--split_antimeridian", action="store_true", default=False
    )
    parser.add_argument("-aggregate", "--aggregate", action="store_true")
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string options for dggrid2geo/dggridgen",
    )
    args = parser.parse_args()

    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    dggrid_instance = create_dggrid_instance()
    result = dggridexpand(
        dggrid_instance=dggrid_instance,
        dggs_type=args.dggs_type,
        input_data=args.input,
        resolution=args.resolution,
        target_resolution=args.target_resolution,
        dggrid_id=args.cellid,
        output_format=args.output_format,
        split_antimeridian=args.split_antimeridian,
        aggregate=args.aggregate,
        options=options,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    dggridcompact_cli()

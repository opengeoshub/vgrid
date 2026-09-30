"""
DGGS vector data to Cloud Optimized GeoParquet.

With -agg_col, coarser levels are parent cells. A resolution-5 row is the
parent of resolution-6 children, with -agg_col summed across those children,
down to -min_res. The file resolution is stored last, unchanged. Without
-agg_col the input is written as one level and -min_res is ignored. Levels
follow

    min(max_res, max(min_res, floor((zoom_level - 2) * 0.8)))

and are written coarse-to-fine so @cogp/reader can fetch only the prefix for
the current zoom. Without -agg_col every input column is copied through and
no resolution column is added.

Key Functions:
    s22cogp: Convert an S2 vector layer into Cloud Optimized GeoParquet
    s22cogp_cli: Command-line interface
"""

import argparse
import json
import math
import os
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from tqdm import tqdm

from vgrid.conversion.dggs2cogp.common import take_rows
from vgrid.conversion.dggs2geo.s22geo import s22geo
from vgrid.dggs import s2
from vgrid.utils.constants import DGGS_TYPES, FIX_ANTIMERIDIAN_CHOICES, ROW_GROUP_ROWS
from vgrid.utils.geometry import degrees_per_pixel
from vgrid.utils.io import (
    add_verbose_argument,
    geoparquet_geometry_columns,
    process_input_data_cogp,
    read_geoparquet_metadata,
    require_column,
)

min_res = DGGS_TYPES["s2"]["min_res"]
max_res = DGGS_TYPES["s2"]["max_res"]
DEFAULT_MIN_RES = min_res
DEFAULT_S2_COL = "s2"
DEFAULT_FIX_ANTIMERIDIAN = "shift"
RES_COLUMN = "res"


def zoom_for_s2_resolution(s2_resolution):
    """Coarsest zoom whose floor((zoom - 2) * 0.8) result is this S2 resolution."""
    return 2.0 + s2_resolution / 0.8


def level_resolution(s2_resolution):
    return degrees_per_pixel(zoom_for_s2_resolution(s2_resolution))


def s2_cell_resolution(token):
    return s2.CellId.from_token(str(token)).level()


def s2_cell_parent(token, resolution):
    return s2.CellId.from_token(str(token)).parent(resolution).to_token()


def detect_source_resolution(cells):
    for cell in cells:
        if cell:
            return s2_cell_resolution(cell)
    raise ValueError("Could not detect an S2 resolution from the file.")


def sum_parents(cells, values, resolution, verbose=True):
    """Sum child values into each parent cell at ``resolution``."""
    totals = {}
    pairs = zip(cells, values)
    for cell, value in tqdm(
        pairs,
        total=len(cells),
        desc=f"Aggregating S2 res {resolution}",
        unit=" cells",
        disable=not verbose,
    ):
        if not cell:
            continue
        parent = s2_cell_parent(cell, resolution)
        totals[parent] = totals.get(parent, 0.0) + value
    return totals


def aggregated_table(
    schema,
    totals,
    resolution,
    s2_col,
    agg_col,
    geometry_col,
    bbox_col,
    fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN,
    verbose=True,
):
    parents = sorted(totals)
    geometries = []
    bboxes = []
    for cell in tqdm(
        parents,
        desc=f"Building S2 res {resolution}",
        unit=" cells",
        disable=not verbose,
    ):
        geom = s22geo(cell, fix_antimeridian=fix_antimeridian)
        if geom is None or getattr(geom, "is_empty", True):
            raise ValueError(f"Could not build a polygon for S2 cell {cell}.")
        geometries.append(geom.wkb)
        bboxes.append(geom.bounds)
    columns = {
        s2_col: pa.array(parents, type=schema.field(s2_col).type),
        RES_COLUMN: pa.array(
            [resolution] * len(parents), type=schema.field(RES_COLUMN).type
        ),
        agg_col: pa.array(
            [totals[cell] for cell in parents], type=schema.field(agg_col).type
        ),
        geometry_col: pa.array(geometries, type=schema.field(geometry_col).type),
    }
    if bbox_col:
        columns[bbox_col] = pa.array(
            [
                {"xmin": bbox[0], "ymin": bbox[1], "xmax": bbox[2], "ymax": bbox[3]}
                for bbox in bboxes
            ],
            type=schema.field(bbox_col).type,
        )
    return pa.table(columns, schema=schema)


def build_lod(level_tables):
    levels = []
    row_group = -1
    for resolution, level_table in level_tables:
        count = level_table.num_rows
        if count:
            row_group += math.ceil(count / ROW_GROUP_ROWS)
        if row_group < 0:
            raise ValueError("The coarsest S2 level has no cells to write.")
        levels.append(
            {
                "row_group_end": row_group,
                "resolution": level_resolution(resolution),
            }
        )
    return levels


def with_lod_metadata(schema, levels, bbox=None, geometry_col=None):
    metadata = dict(schema.metadata or {})
    geo = read_geoparquet_metadata(schema)
    geo["lod"] = {"levels": levels}
    if bbox and geometry_col and geometry_col in geo.get("columns", {}):
        geo["columns"][geometry_col]["bbox"] = [float(value) for value in bbox]
    metadata[b"geo"] = json.dumps(geo, separators=(",", ":")).encode("utf-8")
    return schema.with_metadata(metadata)


def table_extent(level_tables, bbox_col):
    if not bbox_col:
        return None
    xmin = ymin = math.inf
    xmax = ymax = -math.inf
    for _, level_table in level_tables:
        column = level_table.column(bbox_col)
        xmin = min(xmin, pc.min(pc.struct_field(column, "xmin")).as_py())
        ymin = min(ymin, pc.min(pc.struct_field(column, "ymin")).as_py())
        xmax = max(xmax, pc.max(pc.struct_field(column, "xmax")).as_py())
        ymax = max(ymax, pc.max(pc.struct_field(column, "ymax")).as_py())
    if math.isinf(xmin):
        return None
    return [xmin, ymin, xmax, ymax]


def write_cogp(level_tables, output_path, bbox_col, geometry_col):
    levels = build_lod(level_tables)
    schema = with_lod_metadata(
        level_tables[0][1].schema,
        levels,
        table_extent(level_tables, bbox_col),
        geometry_col,
    )
    if os.path.exists(output_path):
        os.remove(output_path)
    writer = pq.ParquetWriter(
        output_path,
        schema,
        compression="zstd",
        compression_level=9,
        write_statistics=True,
        write_page_index=True,
        version="2.6",
    )
    try:
        for level_index, (resolution, level_table) in enumerate(level_tables):
            count = level_table.num_rows
            for start in range(0, count, ROW_GROUP_ROWS):
                writer.write_table(level_table.slice(start, ROW_GROUP_ROWS))
            print(
                f"level {level_index} (s2 resolution {resolution}, "
                f"resolution={level_resolution(resolution):.6f} CRS units): "
                f"{count:,} features"
            )
    finally:
        writer.close()
    return levels


def _bbox_array(values):
    bbox_type = pa.struct(
        [
            pa.field("xmin", pa.float64()),
            pa.field("ymin", pa.float64()),
            pa.field("xmax", pa.float64()),
            pa.field("ymax", pa.float64()),
        ]
    )
    rows = []
    for value in values:
        if hasattr(value, "as_py"):
            value = value.as_py()
        if isinstance(value, dict) and {"xmin", "ymin", "xmax", "ymax"} <= value.keys():
            rows.append(
                {
                    "xmin": value["xmin"],
                    "ymin": value["ymin"],
                    "xmax": value["xmax"],
                    "ymax": value["ymax"],
                }
            )
        else:
            rows.append(None)
    return pa.array(rows, type=bbox_type)


def table_from_geodataframe(gdf):
    """PyArrow table with WKB geometry and GeoParquet metadata."""
    geom_name = gdf.geometry.name
    frame = gdf.drop(columns=geom_name)
    names = []
    arrays = []
    bbox_col = None
    for name in frame.columns:
        names.append(name)
        if name == "geometry_bbox":
            bbox_array = _bbox_array(frame[name].tolist())
            if bbox_array.null_count < len(bbox_array):
                arrays.append(bbox_array)
                bbox_col = name
                continue
        arrays.append(pa.array(frame[name].tolist()))
    wkb = [
        None if geom is None or geom.is_empty else geom.wkb for geom in gdf.geometry
    ]
    names.append(geom_name)
    arrays.append(pa.array(wkb, type=pa.binary()))
    table = pa.Table.from_arrays(arrays, names=names)
    column = {"encoding": "WKB"}
    if bbox_col:
        column["covering"] = {
            "bbox": {
                "xmin": [bbox_col, "xmin"],
                "ymin": [bbox_col, "ymin"],
                "xmax": [bbox_col, "xmax"],
                "ymax": [bbox_col, "ymax"],
            }
        }
    geo = {
        "version": "1.1.0",
        "primary_column": geom_name,
        "columns": {geom_name: column},
    }
    metadata = {b"geo": json.dumps(geo, separators=(",", ":")).encode("utf-8")}
    return table.replace_schema_metadata(metadata)


def input_table(input_path, s2_col, agg_col):
    """Source table for one conversion.

    Without aggregation, Parquet is read directly so every input column and
    its values are kept. Other formats are loaded as a GeoDataFrame first.
    """
    parquet = isinstance(input_path, str) and input_path.lower().endswith(
        (".parquet", ".geoparquet")
    )
    if agg_col is None and parquet:
        table = pq.read_table(input_path)
        require_column(table.column_names, s2_col, "Input")
        return table
    gdf = process_input_data_cogp(input_path, id_col=s2_col, agg_col=agg_col)
    return table_from_geodataframe(gdf)


def output_schema(table, s2_col, agg_col, geometry_col, bbox_col):
    if agg_col:
        names = [s2_col, agg_col]
        if geometry_col not in names:
            names.append(geometry_col)
        if bbox_col and bbox_col not in names:
            names.append(bbox_col)
        for name in names:
            require_column(table.schema.names, name, "Output")
        fields = [table.schema.field(name) for name in names if name != RES_COLUMN]
        s2_index = next(
            index for index, field in enumerate(fields) if field.name == s2_col
        )
        fields.insert(s2_index + 1, pa.field(RES_COLUMN, pa.int8(), nullable=False))
    else:
        fields = list(table.schema)
    return pa.schema(fields, metadata=table.schema.metadata)


def s22cogp(
    input_path,
    output_path,
    agg_col=None,
    s2_col=DEFAULT_S2_COL,
    min_res=DEFAULT_MIN_RES,
    fix_antimeridian=DEFAULT_FIX_ANTIMERIDIAN,
    verbose=True,
):
    """Convert an S2 vector layer into Cloud Optimized GeoParquet.

    Parameters
    ----------
    input_path : str or GeoDataFrame
        GeoParquet, GeoJSON, Shapefile, GeoPackage, or an in-memory layer.
    output_path : str
        Output ``.cogp.parquet`` file.
    agg_col : str, optional
        Numeric column summed into each parent cell. When omitted, only the
        input resolution is written and every input column is kept, with no
        extra resolution column.
    s2_col : str, default "s2"
        S2 cell id column.
    min_res : int, default S2 min_res
        Coarsest parent resolution to build. Ignored when ``agg_col`` is omitted.
    fix_antimeridian : str, default "shift"
        Antimeridian fixing for parent cells: shift, shift_balanced, shift_west,
        shift_east, split, none. Source rows keep the geometry from the input file.
    verbose : bool, default True
        Show progress bars.
    """
    if not agg_col:
        agg_col = None
    if agg_col and s2_col == agg_col:
        raise ValueError("-s2_col and -agg_col must be different columns.")
    if agg_col and (s2_col == RES_COLUMN or agg_col == RES_COLUMN):
        raise ValueError(f"'{RES_COLUMN}' is reserved for the S2 resolution column.")
    if fix_antimeridian == "none":
        fix_antimeridian = None
    elif fix_antimeridian is not None and fix_antimeridian not in FIX_ANTIMERIDIAN_CHOICES:
        raise ValueError(
            "-fix_antimeridian must be one of "
            + ", ".join(FIX_ANTIMERIDIAN_CHOICES)
            + f", got {fix_antimeridian}."
        )

    table = input_table(input_path, s2_col, agg_col)

    cells = [
        None if value is None else str(value) for value in table.column(s2_col).to_pylist()
    ]
    source_resolution = detect_source_resolution(cells)
    if not DEFAULT_MIN_RES <= source_resolution <= max_res:
        raise ValueError(
            f"S2 resolution must be between {DEFAULT_MIN_RES} and "
            f"{max_res}, got {source_resolution}."
        )
    if agg_col and not DEFAULT_MIN_RES <= min_res <= source_resolution:
        raise ValueError(
            f"-min_res must be between {DEFAULT_MIN_RES} and the file "
            f"resolution {source_resolution}, got {min_res}."
        )

    geo = read_geoparquet_metadata(table.schema)
    geometry_col, bbox_col = geoparquet_geometry_columns(geo)
    require_column(table.schema.names, geometry_col, "Geometry")
    schema = output_schema(table, s2_col, agg_col, geometry_col, bbox_col)

    level_tables = []
    if agg_col:
        values = [
            0.0 if value is None else float(value)
            for value in table.column(agg_col).to_pylist()
        ]
        print(
            f"Aggregating {len(cells):,} on column '{agg_col}' at S2 resolution "
            f"{source_resolution} into resolutions {min_res}..{source_resolution}"
        )
        for resolution in range(min_res, source_resolution):
            totals = sum_parents(cells, values, resolution, verbose=verbose)
            level_tables.append(
                (
                    resolution,
                    aggregated_table(
                        schema,
                        totals,
                        resolution,
                        s2_col,
                        agg_col,
                        geometry_col,
                        bbox_col,
                        fix_antimeridian=fix_antimeridian,
                        verbose=verbose,
                    ),
                )
            )
    else:
        print(
            f"Writing {len(cells):,} S2 cells at resolution {source_resolution}"
        )
    order = sorted(range(len(cells)), key=lambda index: cells[index] or "")
    row_index = pa.array(order, type=pa.int64())
    if agg_col:
        source_names = [name for name in schema.names if name != RES_COLUMN]
        source = take_rows(table.select(source_names), row_index)
        source = source.append_column(
            pa.field(RES_COLUMN, pa.int8(), nullable=False),
            pa.array([source_resolution] * source.num_rows, type=pa.int8()),
        ).select(schema.names)
    else:
        source = take_rows(table.select(schema.names), row_index)
    level_tables.append((source_resolution, source))
    write_cogp(level_tables, output_path, bbox_col, geometry_col)
    print(f"Wrote {output_path}")


def s22cogp_cli():
    parser = argparse.ArgumentParser(
        description=(
            "Convert an S2 vector layer to Cloud Optimized GeoParquet. "
            "With -agg_col, each coarser level is parent cells summed from their children."
        )
    )
    parser.add_argument(
        "input",
        help="Input file, for example s2_16.parquet, s2_16.geojson, or s2_16.gpkg.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output .cogp.parquet file. Default: <input>.cogp.parquet.",
    )
    parser.add_argument(
        "-s2_col",
        default=DEFAULT_S2_COL,
        help=f"S2 cell column. Default: {DEFAULT_S2_COL}.",
    )
    parser.add_argument(
        "-agg_col",
        default=None,
        help=(
            "Numeric column to sum into each parent cell. "
            "Omit to write the input resolution only, keeping every input column "
            "and not adding a resolution column."
        ),
    )
    parser.add_argument(
        "-min_res",
        type=int,
        default=DEFAULT_MIN_RES,
        help=(
            "Coarsest parent resolution to build when -agg_col is set. "
            "Resolution 5 sums resolution-6 children, and so on down to this value. "
            f"Ignored without -agg_col. Default: {DEFAULT_MIN_RES}."
        ),
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=DEFAULT_FIX_ANTIMERIDIAN,
        help=(
            "Antimeridian fixing method: shift, shift_balanced, shift_west, "
            f"shift_east, split, none. Default: {DEFAULT_FIX_ANTIMERIDIAN}."
        ),
    )
    add_verbose_argument(parser)
    args = parser.parse_args()
    output_path = args.output
    if not output_path:
        base = os.path.splitext(args.input)[0]
        output_path = f"{base}.cogp.parquet"
    try:
        s22cogp(
            args.input,
            output_path,
            args.agg_col,
            args.s2_col,
            args.min_res,
            fix_antimeridian=args.fix_antimeridian,
            verbose=args.verbose,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    s22cogp_cli()

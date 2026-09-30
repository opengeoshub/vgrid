"""
Shared engine for the ``<dggs>2cogp`` modules.

A ``CogpSpec`` describes one DGGS: how to read a cell's resolution, how to map
cells to parents and parents to polygons (both batched), and which resolution
Vgrid Viz shows at a given zoom. ``dggs2cogp`` does the rest the same way as
``h32cogp``: with an aggregation column, coarser levels are parent cells with
the column summed from their children, down to ``min_res``; the file
resolution is written last. Levels are written coarse-to-fine and each level's
LOD ``resolution`` is the degrees-per-pixel of the coarsest zoom at which
Vgrid Viz switches to that DGGS resolution, so @cogp/reader can fetch only
the prefix for the current zoom.
"""

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from typing import Callable, List, Optional

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from tqdm import tqdm

from vgrid.utils.constants import FIX_ANTIMERIDIAN_CHOICES, ROW_GROUP_ROWS
from vgrid.utils.geometry import degrees_per_pixel
from vgrid.utils.io import (
    add_verbose_argument,
    geoparquet_geometry_columns,
    process_input_data_cogp,
    read_geoparquet_metadata,
    require_column,
)

RES_COLUMN = "res"
# Same zoom <-> scale denominator relation as the Vgrid Viz layers.
_ZOOM_SCALE_OFFSET = 29.1402
_MAX_ZOOM = 35.0
_ZOOM_STEP = 0.01


def scale_for_zoom(zoom):
    return 2.0 ** (_ZOOM_SCALE_OFFSET - zoom)


_WIDE_TYPES = {pa.binary(): pa.large_binary(), pa.string(): pa.large_string()}
_NARROW_TYPES = {wide: narrow for narrow, wide in _WIDE_TYPES.items()}


def _map_field_types(schema, mapping):
    return pa.schema(
        [field.with_type(mapping.get(field.type, field.type)) for field in schema],
        metadata=schema.metadata,
    )


def take_rows(table, indices):
    """``table.take(indices)`` that works past 2 GB of binary/string data.

    ``take`` concatenates chunks, which overflows 32-bit offsets on large
    geometry columns, so values are widened to 64-bit offsets first.
    """
    return table.cast(_map_field_types(table.schema, _WIDE_TYPES)).take(indices)


def narrow_schema(schema):
    """Undo ``take_rows`` widening so the written Parquet schema is unchanged."""
    return _map_field_types(schema, _NARROW_TYPES)


@dataclass
class CogpSpec:
    label: str
    min_res: int
    max_res: int
    cell_resolution: Optional[Callable[[str], int]]
    parents: Callable[[List[str], int], List[str]]
    geometries: Callable[[List[str], int], list]
    resolution_for_zoom: Callable[[float], int]
    valid_resolutions: Optional[List[int]] = None

    def levels_between(self, low, high):
        """Valid resolutions in ``[low, high)``."""
        if self.valid_resolutions is None:
            return list(range(low, high))
        return [r for r in self.valid_resolutions if low <= r < high]

    def is_valid(self, resolution):
        if self.valid_resolutions is not None:
            return resolution in self.valid_resolutions
        return self.min_res <= resolution <= self.max_res


def per_cell_parents(parent_at):
    """Lift ``parent_at(cell, resolution)`` to the batched ``parents`` signature."""
    return lambda cells, resolution: [parent_at(c, resolution) for c in cells]


def per_cell_geometries(geo_fn):
    """Lift ``geo_fn(cell, resolution)`` to the batched ``geometries`` signature."""
    return lambda cells, resolution: [geo_fn(c, resolution) for c in cells]


def zoom_for_resolution(spec, resolution):
    """Coarsest zoom at which ``spec.resolution_for_zoom`` reaches ``resolution``."""
    return _zoom_for_resolution(spec.resolution_for_zoom, resolution)


@lru_cache(maxsize=None)
def _zoom_for_resolution(resolution_for_zoom, resolution):
    steps = int(_MAX_ZOOM / _ZOOM_STEP)
    for step in range(steps + 1):
        zoom = step * _ZOOM_STEP
        if resolution_for_zoom(zoom) >= resolution:
            return zoom
    return _MAX_ZOOM


def level_resolution(spec, resolution):
    return degrees_per_pixel(zoom_for_resolution(spec, resolution))


def sum_parents(spec, cells, values, resolution, verbose=True):
    """Sum child values into each parent cell at ``resolution``."""
    present = [(c, v) for c, v in zip(cells, values) if c]
    parents = spec.parents([c for c, _ in present], resolution)
    totals = {}
    for parent, (_, value) in tqdm(
        zip(parents, present),
        total=len(present),
        desc=f"Aggregating {spec.label} res {resolution}",
        unit=" cells",
        disable=not verbose,
    ):
        totals[parent] = totals.get(parent, 0.0) + value
    return totals


def aggregated_table(
    spec,
    schema,
    totals,
    resolution,
    id_col,
    agg_col,
    geometry_col,
    bbox_col,
    verbose=True,
):
    parents = sorted(totals)
    if verbose:
        print(f"Building {len(parents):,} {spec.label} res {resolution} polygons")
    geoms = spec.geometries(parents, resolution)
    geometries = []
    bboxes = []
    for cell, geom in zip(parents, geoms):
        if geom is None or getattr(geom, "is_empty", True):
            raise ValueError(f"Could not build a polygon for {spec.label} cell {cell}.")
        geometries.append(geom.wkb)
        bboxes.append(geom.bounds)
    columns = {
        id_col: pa.array(parents, type=schema.field(id_col).type),
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


def build_lod(spec, level_tables):
    levels = []
    row_group = -1
    for resolution, level_table in level_tables:
        count = level_table.num_rows
        if count:
            row_group += math.ceil(count / ROW_GROUP_ROWS)
        if row_group < 0:
            raise ValueError(f"The coarsest {spec.label} level has no cells to write.")
        levels.append(
            {
                "row_group_end": row_group,
                "resolution": level_resolution(spec, resolution),
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


def write_cogp(spec, level_tables, output_path, bbox_col, geometry_col):
    levels = build_lod(spec, level_tables)
    schema = with_lod_metadata(
        narrow_schema(level_tables[0][1].schema),
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
                writer.write_table(
                    level_table.slice(start, ROW_GROUP_ROWS).cast(schema)
                )
            print(
                f"level {level_index} ({spec.label} resolution {resolution}, "
                f"resolution={level_resolution(spec, resolution):.6f} CRS units): "
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
    wkb = [None if geom is None or geom.is_empty else geom.wkb for geom in gdf.geometry]
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


def input_table(input_path, id_col, agg_col):
    """Source table. Without aggregation, Parquet is read directly so every
    input column and its values are kept."""
    parquet = isinstance(input_path, str) and input_path.lower().endswith(
        (".parquet", ".geoparquet")
    )
    if agg_col is None and parquet:
        table = pq.read_table(input_path)
        require_column(table.column_names, id_col, "Input")
        return table
    gdf = process_input_data_cogp(input_path, id_col=id_col, agg_col=agg_col)
    return table_from_geodataframe(gdf)


def output_schema(table, id_col, agg_col, geometry_col, bbox_col):
    if agg_col:
        names = [id_col, agg_col]
        if geometry_col not in names:
            names.append(geometry_col)
        if bbox_col and bbox_col not in names:
            names.append(bbox_col)
        for name in names:
            require_column(table.schema.names, name, "Output")
        fields = [table.schema.field(name) for name in names if name != RES_COLUMN]
        id_index = next(
            index for index, field in enumerate(fields) if field.name == id_col
        )
        fields.insert(id_index + 1, pa.field(RES_COLUMN, pa.int8(), nullable=False))
    else:
        fields = list(table.schema)
    return pa.schema(fields, metadata=table.schema.metadata)


def normalize_fix_antimeridian(fix_antimeridian):
    if fix_antimeridian == "none":
        return None
    if fix_antimeridian is not None and fix_antimeridian not in FIX_ANTIMERIDIAN_CHOICES:
        raise ValueError(
            "-fix_antimeridian must be one of "
            + ", ".join(FIX_ANTIMERIDIAN_CHOICES)
            + f", got {fix_antimeridian}."
        )
    return fix_antimeridian


def dggs2cogp(
    spec,
    input_path,
    output_path,
    agg_col=None,
    id_col=None,
    min_res=None,
    source_resolution=None,
    verbose=True,
):
    """Convert a DGGS vector layer described by ``spec`` into Cloud Optimized GeoParquet.

    ``source_resolution`` is required when ``spec.cell_resolution`` is None
    (DGGRID SEQNUMs); otherwise it is read from the first cell.
    """
    label = spec.label
    if min_res is None:
        min_res = spec.min_res
    if not agg_col:
        agg_col = None
    if agg_col and id_col == agg_col:
        raise ValueError("The cell column and -agg_col must be different columns.")
    if agg_col and (id_col == RES_COLUMN or agg_col == RES_COLUMN):
        raise ValueError(f"'{RES_COLUMN}' is reserved for the {label} resolution column.")

    table = input_table(input_path, id_col, agg_col)
    cells = [
        None if value is None else str(value)
        for value in table.column(id_col).to_pylist()
    ]
    if source_resolution is None:
        if spec.cell_resolution is None:
            raise ValueError(f"{label} needs the input resolution.")
        first = next((c for c in cells if c), None)
        if first is None:
            raise ValueError(f"Could not detect a {label} resolution from the file.")
        source_resolution = spec.cell_resolution(first)
    if not spec.is_valid(source_resolution):
        raise ValueError(
            f"{label} resolution must be between {spec.min_res} and "
            f"{spec.max_res}, got {source_resolution}."
        )
    if agg_col and not (spec.min_res <= min_res <= source_resolution):
        raise ValueError(
            f"-min_res must be between {spec.min_res} and the file "
            f"resolution {source_resolution}, got {min_res}."
        )

    geo = read_geoparquet_metadata(table.schema)
    geometry_col, bbox_col = geoparquet_geometry_columns(geo)
    require_column(table.schema.names, geometry_col, "Geometry")
    schema = output_schema(table, id_col, agg_col, geometry_col, bbox_col)

    level_tables = []
    if agg_col:
        values = [
            0.0 if value is None else float(value)
            for value in table.column(agg_col).to_pylist()
        ]
        levels = spec.levels_between(min_res, source_resolution)
        print(
            f"Aggregating {len(cells):,} on column '{agg_col}' at {label} resolution "
            f"{source_resolution} into resolutions {levels + [source_resolution]}"
        )
        for resolution in levels:
            totals = sum_parents(spec, cells, values, resolution, verbose=verbose)
            level_tables.append(
                (
                    resolution,
                    aggregated_table(
                        spec,
                        schema,
                        totals,
                        resolution,
                        id_col,
                        agg_col,
                        geometry_col,
                        bbox_col,
                        verbose=verbose,
                    ),
                )
            )
    else:
        print(f"Writing {len(cells):,} {label} cells at resolution {source_resolution}")
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
    write_cogp(spec, level_tables, output_path, bbox_col, geometry_col)
    print(f"Wrote {output_path}")


def cogp_arg_parser(label, col_flag, default_col, default_min_res, example):
    """Parser with the options every ``<dggs>2cogp`` CLI shares."""
    parser = argparse.ArgumentParser(
        description=(
            f"Convert a {label} vector layer to Cloud Optimized GeoParquet. "
            "With -agg_col, each coarser level is parent cells summed from their children."
        )
    )
    parser.add_argument(
        "input",
        help=f"Input file, for example {example}.parquet, {example}.geojson, or {example}.gpkg.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output .cogp.parquet file. Default: <input>.cogp.parquet.",
    )
    parser.add_argument(
        f"-{col_flag}",
        dest="id_col",
        default=default_col,
        help=f"{label} cell column. Default: {default_col}.",
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
        default=default_min_res,
        help=(
            "Coarsest parent resolution to build when -agg_col is set. "
            f"Ignored without -agg_col. Default: {default_min_res}."
        ),
    )
    add_verbose_argument(parser)
    return parser


def add_fix_antimeridian_argument(parser, default):
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=default,
        help=(
            "Antimeridian fixing method for parent cells: shift, shift_balanced, "
            f"shift_west, shift_east, split, none. Default: {default}."
        ),
    )


def add_split_antimeridian_argument(parser):
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Split parent cells at the antimeridian.",
    )


def default_output_path(args):
    if args.output:
        return args.output
    return f"{os.path.splitext(args.input)[0]}.cogp.parquet"


def run_cli(convert):
    """Call ``convert()`` and exit with status 1 on error, like ``h32cogp_cli``."""
    try:
        convert()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

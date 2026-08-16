"""
A5 Segments Analysis Module

Evaluate how A5 cell polygon geometry/area changes with the ``segments`` option
passed to ``a52geo`` / ``a5.cell_to_boundary``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import a5
import geopandas as gpd
from a5.core.cell_info import get_num_cells
from tqdm import tqdm

from vgrid.conversion.dggs2geo.a52geo import a52geo
from vgrid.utils.constants import AUTHALIC_AREA
from vgrid.utils.geometry import geod

DEFAULT_A5_ID = "6361180000000000"
DEFAULT_SEGMENTS = (1, 10, 100, 1000, 10000, 100000, 1000000)


def a5_segments(
    a5_id: str = DEFAULT_A5_ID,
    segments=None,
    output: str | None = None,
    split_antimeridian: bool = False,
) -> gpd.GeoDataFrame:
    """
    Build A5 polygons for a single cell across a sequence of ``segments`` values.

    Parameters
    ----------
    a5_id : str, default ``'6361180000000000'``
        A5 cell ID (hex string).
    segments : sequence of int, optional
        Segment counts to test. Defaults to ``1, 10, 100, ..., 1000000``.
    output : str, optional
        Output GeoParquet path. Defaults to ``a5_segments_{a5_id}.parquet``.
    split_antimeridian : bool, default False
        Passed through to ``a52geo``.

    Returns
    -------
    geopandas.GeoDataFrame
        Columns: ``a5_id``, ``resolution``, ``segments``, ``cell_area``,
        ``norm_area``, ``area_error``, ``geometry``.

        ``area_error`` is ``100 * abs(norm_area - 1)`` (percent).
    """
    if segments is None:
        segments = DEFAULT_SEGMENTS

    resolution = a5.get_resolution(a5.hex_to_u64(a5_id))
    mean_area = AUTHALIC_AREA / get_num_cells(resolution)

    rows = []
    for n_segments in tqdm(list(segments), desc="A5 segments", unit=" value"):
        cell_polygon = a52geo(
            a5_id,
            options={"segments": int(n_segments)},
            split_antimeridian=split_antimeridian,
        )
        if cell_polygon is None or getattr(cell_polygon, "is_empty", False):
            continue

        cell_area_perimeter = geod.geometry_area_perimeter(cell_polygon)
        cell_area = abs(cell_area_perimeter[0])
        norm_area = cell_area / mean_area
        rows.append(
            {
                "a5_id": a5_id,
                "resolution": resolution,
                "segments": int(n_segments),
                "cell_area": cell_area,
                "norm_area": norm_area,
                "area_error": 100 * abs(norm_area - 1),
                "geometry": cell_polygon,
            }
        )

    if not rows:
        raise ValueError(f"No valid A5 polygons produced for a5_id={a5_id!r}")

    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    if output is None:
        output = f"a5_segments.parquet"
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(output_path, index=False)
    print(f"Saved {len(gdf)} rows to {output_path.resolve()}")
    return gdf


def a5_segments_cli():
    """Command-line interface for A5 segments analysis."""
    parser = argparse.ArgumentParser(
        description="Evaluate A5 cell_area across a52geo segments values"
    )
    parser.add_argument(
        "-a5",
        "--a5_id",
        type=str,
        default=DEFAULT_A5_ID,
        help=f"A5 cell ID (default: {DEFAULT_A5_ID})",
    )
    parser.add_argument(
        "-s",
        "--segments",
        type=int,
        nargs="+",
        default=list(DEFAULT_SEGMENTS),
        help="Segment counts to evaluate (default: 1 10 100 ... 1000000)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output GeoParquet path",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable antimeridian splitting",
    )
    args = parser.parse_args()
    a5_segments(
        a5_id=args.a5_id,
        segments=args.segments,
        output=args.output,
        split_antimeridian=args.split_antimeridian,
    )


if __name__ == "__main__":
    a5_segments_cli()

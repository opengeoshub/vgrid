"""
DGGS resampling (area-weighted or nearest-neighbour attribute transfer).

Port of the QGIS vgridtools ``dggsresample`` workflow: pick a target resolution
(optionally by matching mean cell area), build a target DGGS grid over the
source footprint, then transfer a numeric field using either area-weighted overlap or
nearest-neighbour assignment from the nearest source cell centroid.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import sys
from numbers import Number
from typing import Any, Optional, Union

import geopandas as gpd
import h3
import a5
import pandas as pd
from tqdm import tqdm

from vgrid.dggs import s2, olc, mercantile
from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.generator.a5grid import a5_grid
from vgrid.generator.geohashgrid import geohash_grid_within_bbox
from vgrid.generator.h3grid import h3_grid_within_bbox
from vgrid.generator.olcgrid import olc_grid_within_bbox
from vgrid.generator.quadkeygrid import quadkey_grid
from vgrid.generator.qtmgrid import qtm_grid_within_bbox
from vgrid.generator.rhealpixgrid import rhealpix_grid_within_bbox
from vgrid.generator.healpixgrid import healpix_grid_within_bbox
from vgrid.generator.s2grid import s2_grid
from vgrid.generator.tilecodegrid import tilecode_grid
from vgrid.stats.a5stats import a5_metrics
from vgrid.stats.geohashstats import geohash_metrics
from vgrid.stats.isea4tstats import isea4t_metrics
from vgrid.stats.olcstats import olc_metrics
from vgrid.stats.quadkeystats import quadkey_metrics
from vgrid.stats.qtmstats import qtm_metrics
from vgrid.stats.rhealpixstats import rhealpix_metrics
from vgrid.stats.s2stats import s2_metrics
from vgrid.stats.tilecodestats import tilecode_metrics
from vgrid.stats.dggridstats import dggridstats
from vgrid.utils.constants import (
    DGGAL_TYPES,
    DGGRID_TYPES,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)
from shapely import make_valid
from shapely.geometry.base import BaseGeometry

from vgrid.utils.geometry import _reproject_for_metric
from vgrid.utils.io import (
    convert_to_output_format,
    create_dggrid_instance,
    process_input_data_resample,
    validate_bbox,
    validate_dggal_type,
    validate_dggrid_type,
    olc_resolutions,
)
from vgrid.stats.dggalstats import dggal_metrics
from vgrid.generator.dggalgen import dggalgen
from vgrid.generator.dggridgen import generate_grid as dggrid_generate_polygons
from vgrid.conversion.dggscompact.dggridcompact import _resolve_cell_geometry

import dggal

E = WGS84_ELLIPSOID

if platform.system() == "Windows":
    from vgrid.generator.isea4tgrid import isea4t_grid_within_bbox


def _dggal_short_type(dggs_key: str) -> Optional[str]:
    """
    Map ``dggal_gnosis`` / ``dggal_rtea3h`` / … to keys in ``DGGAL_TYPES``.

    Bare names in ``DGGAL_TYPES`` are accepted except ``rhealpix``, which stays
    mapped to vgrid's native rHEALPix grid — use ``dggal_rhealpix`` for DGGAL.
    """
    k = dggs_key.strip().lower()
    if k.startswith("dggal_"):
        short = k[len("dggal_") :]
        return short if short in DGGAL_TYPES else None
    if k in DGGAL_TYPES and k != "rhealpix":
        return k
    return None


def _dggrid_short_type(dggs_key: str) -> Optional[str]:
    """
    Map ``dggrid_ISEA4T`` / ``dggrid_isea7h`` / … to keys in ``DGGRID_TYPES``.

    Bare identifiers must match ``DGGRID_TYPES`` exactly (e.g. ``ISEA7H``), so
    they do not collide with native lowercase keys such as ``isea4t``.
    """
    raw = dggs_key.strip()
    low = raw.lower()
    if low.startswith("dggrid_"):
        short = raw[len("dggrid_") :].strip().upper()
        return short if short in DGGRID_TYPES else None
    if raw in DGGRID_TYPES:
        return raw
    return None


def get_nearest_resolution(
    source_gdf: gpd.GeoDataFrame,
    from_dggs: str,
    to_dggs: str,
    from_col: Optional[str] = None,
) -> int:
    """
    Match mean cell area of the first source cell to the closest ``to_dggs``
    resolution (same search as the QGIS plugin).
    """
    if from_col is None:
        from_col = from_dggs

    if source_gdf.empty or from_col not in source_gdf.columns:
        raise ValueError(f"No features or missing <{from_col}> column.")

    from_dggs_id = source_gdf.iloc[0][from_col]
    if from_dggs_id is None or (
        isinstance(from_dggs_id, float) and pd.isna(from_dggs_id)
    ):
        raise ValueError(f"No valid DGGS IDs found in <{from_col}> column.")

    try:
        if from_dggs == "h3":
            from_resolution = h3.get_resolution(from_dggs_id)
            from_area = h3.average_hexagon_area(from_resolution, unit="m^2")

        elif from_dggs == "s2":
            s2_id = s2.CellId.from_token(from_dggs_id)
            from_resolution = s2_id.level()
            _, _, from_area, _ = s2_metrics(from_resolution)

        elif from_dggs == "a5":
            from_resolution = a5.get_resolution(a5.hex_to_u64(from_dggs_id))
            _, _, from_area, _ = a5_metrics(from_resolution)

        elif from_dggs == "rhealpix":
            rhealpix_uids = (from_dggs_id[0],) + tuple(map(int, from_dggs_id[1:]))
            rhealpix_dggs = RHEALPixDGGS(
                ellipsoid=E, north_square=1, south_square=3, N_side=3
            )
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            from_resolution = rhealpix_cell.resolution
            _, _, from_area, _ = rhealpix_metrics(from_resolution)

        elif from_dggs == "isea4t":
            if platform.system() != "Windows":
                raise ValueError(
                    "isea4t area matching is only available on Windows in this build.",
                )
            from_resolution = len(from_dggs_id) - 2
            _, _, from_area, _ = isea4t_metrics(from_resolution)

        elif from_dggs == "qtm":
            from_resolution = len(from_dggs_id)
            _, _, from_area, _ = qtm_metrics(from_resolution)

        elif from_dggs == "olc":
            coord = olc.decode(from_dggs_id)
            from_resolution = coord.codeLength
            _, _, from_area, _ = olc_metrics(from_resolution)

        elif from_dggs == "geohash":
            from_resolution = len(from_dggs_id)
            _, _, from_area, _ = geohash_metrics(from_resolution)

        elif from_dggs == "tilecode":
            match = re.match(r"z(\d+)x(\d+)y(\d+)", from_dggs_id)
            if not match:
                raise ValueError("tilecode id must look like z{x}x{y}y{y}")
            from_resolution = int(match.group(1))
            _, _, from_area, _ = tilecode_metrics(from_resolution)

        elif from_dggs == "quadkey":
            tile = mercantile.quadkey_to_tile(from_dggs_id)
            from_resolution = tile.z
            _, _, from_area, _ = quadkey_metrics(from_resolution)

        elif (dt := _dggal_short_type(from_dggs)) is not None:
            dt = validate_dggal_type(dt)
            cls_name = DGGAL_TYPES[dt]["class_name"]
            dggrs = getattr(dggal, cls_name)()
            zone = dggrs.getZoneFromTextID(str(from_dggs_id))
            from_resolution = int(dggrs.getZoneLevel(zone))
            _, _, from_area, _ = dggal_metrics(dt, from_resolution)

        elif (dt := _dggrid_short_type(from_dggs)) is not None:
            dt = validate_dggrid_type(dt)
            dggrid_instance = create_dggrid_instance()
            if "resolution" in source_gdf.columns:
                rv = source_gdf.iloc[0]["resolution"]
                if rv is None or (isinstance(rv, float) and pd.isna(rv)):
                    raise ValueError(
                        "Missing or invalid <resolution> for DGGRID source."
                    )
                from_resolution = int(rv)
            else:
                pref = int(DGGRID_TYPES[dt]["default_res"])
                _, from_resolution = _resolve_cell_geometry(
                    dggrid_instance,
                    dt,
                    str(from_dggs_id),
                    preferred_resolution=pref,
                )
                if from_resolution is None:
                    raise ValueError(
                        f"Could not resolve DGGRID cell id {from_dggs_id!r} for type {dt}."
                    )
            stats = dggridstats(dggrid_instance, dt, unit="m")
            row = stats[stats["resolution"] == from_resolution]
            if row.empty:
                raise ValueError(
                    f"No DGGRID stats for {dt} at resolution {from_resolution}."
                )
            from_area = float(row["area_m2"].iloc[0])

        else:
            raise ValueError(f"Unsupported from_dggs type for area match: {from_dggs}")

    except Exception as e:
        raise ValueError(f"Failed to calculate area from {from_dggs}: {e}") from e

    nearest_resolution: int
    min_diff = float("inf")

    try:
        if to_dggs == "h3":
            for res in range(16):
                avg_area = h3.average_hexagon_area(res, unit="m^2")
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "s2":
            for res in range(31):
                _, _, avg_area, _ = s2_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "a5":
            for res in range(30):
                _, _, avg_area, _ = a5_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "rhealpix":
            for res in range(16):
                _, _, avg_area, _ = rhealpix_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "isea4t":
            if platform.system() != "Windows":
                raise ValueError(
                    "isea4t nearest resolution is only available on Windows in this build.",
                )
            for res in range(26):
                _, _, avg_area, _ = isea4t_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "qtm":
            for res in range(1, 25):
                _, _, avg_area, _ = qtm_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "olc":
            for res in olc_resolutions:
                _, _, avg_area, _ = olc_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "geohash":
            for res in range(1, 11):
                _, _, avg_area, _ = geohash_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "tilecode":
            for res in range(30):
                _, _, avg_area, _ = tilecode_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif to_dggs == "quadkey":
            for res in range(30):
                _, _, avg_area, _ = quadkey_metrics(res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif (dt := _dggal_short_type(to_dggs)) is not None:
            dt = validate_dggal_type(dt)
            lo = int(DGGAL_TYPES[dt]["min_res"])
            hi = int(DGGAL_TYPES[dt]["max_res"])
            for res in range(lo, hi + 1):
                _, _, avg_area, _ = dggal_metrics(dt, res)
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        elif (dt := _dggrid_short_type(to_dggs)) is not None:
            dt = validate_dggrid_type(dt)
            dggrid_instance = create_dggrid_instance()
            stats = dggridstats(dggrid_instance, dt, unit="m")
            lo = int(DGGRID_TYPES[dt]["min_res"])
            hi = int(DGGRID_TYPES[dt]["max_res"])
            for _, row in stats.iterrows():
                res = int(row["resolution"])
                if res < lo or res > hi:
                    continue
                avg_area = float(row["area_m2"])
                diff = abs(avg_area - from_area)
                if diff < min_diff:
                    min_diff = diff
                    nearest_resolution = res

        else:
            raise ValueError(f"Unsupported to_dggs type for area match: {to_dggs}")

    except Exception as e:
        raise ValueError(
            f"Failed to calculate nearest resolution for {to_dggs}: {e}"
        ) from e

    return nearest_resolution


def generate_grid(
    source_gdf: gpd.GeoDataFrame,
    to_dggs: str,
    resolution: int,
    *,
    fix_antimeridian: Optional[str] = None,
    split_antimeridian: bool = False,
    aggregate: bool = False,
    dggrid_options: Optional[dict] = None,
    a5_options: Optional[dict] = None,
) -> gpd.GeoDataFrame:
    """
    Build a target DGGS GeoDataFrame over the axis-aligned bounds of ``source_gdf``
    (``source_gdf.total_bounds``), using bbox-scoped grid generators only.
    """
    if source_gdf.empty:
        raise ValueError("No features provided for grid generation.")

    bbox = validate_bbox(list(source_gdf.total_bounds))

    if to_dggs == "h3":
        gdf = h3_grid_within_bbox(resolution, bbox, fix_antimeridian=fix_antimeridian)
    elif to_dggs == "s2":
        gdf = s2_grid(resolution, bbox, fix_antimeridian=fix_antimeridian)
    elif to_dggs == "a5":
        gdf = a5_grid(
            resolution,
            bbox,
            options=a5_options,
            split_antimeridian=split_antimeridian,
        )
    elif to_dggs == "rhealpix":
        gdf = rhealpix_grid_within_bbox(
            resolution, bbox, fix_antimeridian=fix_antimeridian
        )
    elif to_dggs == "healpix":
        gdf = healpix_grid_within_bbox(
            resolution, bbox, fix_antimeridian=fix_antimeridian
        )
    elif to_dggs == "isea4t":
        if platform.system() != "Windows":
            raise ValueError("isea4t grid generation requires Windows in this build.")
        gdf = isea4t_grid_within_bbox(
            resolution, bbox, fix_antimeridian=fix_antimeridian
        )
    elif to_dggs == "qtm":
        gdf = qtm_grid_within_bbox(resolution, bbox)
    elif to_dggs == "olc":
        gdf = olc_grid_within_bbox(resolution, bbox)
    elif to_dggs == "geohash":
        gdf = geohash_grid_within_bbox(resolution, bbox)
    elif to_dggs == "tilecode":
        gdf = tilecode_grid(resolution, bbox)
    elif to_dggs == "quadkey":
        gdf = quadkey_grid(resolution, bbox)
    elif (dt := _dggal_short_type(to_dggs)) is not None:
        dt = validate_dggal_type(dt)
        bbox_t = (float(bbox[0]), float(bbox[1]), float(bbox[2]), float(bbox[3]))
        use_split = split_antimeridian or bool(fix_antimeridian)
        gdf = dggalgen(
            dggs_type=dt,
            resolution=resolution,
            bbox=bbox_t,
            compact=False,
            output_format="gpd",
            split_antimeridian=use_split,
        )
        if gdf is None or getattr(gdf, "empty", True):
            raise ValueError(
                f"No DGGAL cells for type {to_dggs!r} ({dt}), resolution {resolution}, bbox {bbox_t}."
            )
    elif (dt := _dggrid_short_type(to_dggs)) is not None:
        dt = validate_dggrid_type(dt)
        dggrid_instance = create_dggrid_instance()
        use_split = split_antimeridian or bool(fix_antimeridian)
        gdf = dggrid_generate_polygons(
            dggrid_instance,
            dt,
            resolution,
            bbox,
            output_address_type="SEQNUM",
            split_antimeridian=use_split,
            aggregate=aggregate,
            options=dggrid_options,
        )
        if gdf is None or getattr(gdf, "empty", True):
            raise ValueError(
                f"No DGGRID cells for type {to_dggs!r} ({dt}), resolution {resolution}, bbox {bbox}."
            )
    else:
        raise ValueError(f"Unsupported DGGS type: {to_dggs}")

    return gdf


def _ensure_valid_geometry(geom: BaseGeometry) -> BaseGeometry:
    """Repair invalid metric geometries before intersection (e.g. after antimeridian shift)."""
    if geom is None or geom.is_empty:
        return geom
    if geom.is_valid:
        return geom
    fixed = make_valid(geom)
    if fixed.is_empty:
        return fixed
    if fixed.geom_type == "GeometryCollection":
        polys = [
            g
            for g in fixed.geoms
            if g.geom_type in ("Polygon", "MultiPolygon") and not g.is_empty
        ]
        if polys:
            fixed = max(polys, key=lambda g: g.area)
    if not fixed.is_valid:
        fixed = fixed.buffer(0)
    return fixed


def _prepare_metric_geometries(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    out = gdf.copy()
    out.geometry = out.geometry.map(_ensure_valid_geometry)
    return out


def _safe_intersection(geom_a: BaseGeometry, geom_b: BaseGeometry) -> BaseGeometry:
    """Intersection with geometry repair on topological failures."""
    try:
        inter = geom_a.intersection(geom_b)
        if not inter.is_empty:
            return inter
    except Exception:
        pass
    a = _ensure_valid_geometry(geom_a)
    b = _ensure_valid_geometry(geom_b)
    return a.intersection(b)


def _resampling_area_weighted(
    source_gdf: gpd.GeoDataFrame,
    target_gdf: gpd.GeoDataFrame,
    resample_col: str,
) -> gpd.GeoDataFrame:
    if source_gdf.empty or target_gdf.empty:
        return gpd.GeoDataFrame(columns=target_gdf.columns, crs=target_gdf.crs)

    for val in source_gdf[resample_col]:
        if not isinstance(val, Number):
            raise TypeError(
                f"Non-numeric value found in <{resample_col}>. "
                "Resampled field calculation failed."
            )

    source_metric, target_metric = _reproject_for_metric(source_gdf, target_gdf)
    source_metric = _prepare_metric_geometries(source_metric)
    target_metric = _prepare_metric_geometries(target_metric)

    # Spatial index: only evaluate intersections for overlapping pairs (not n × m).
    joined = gpd.sjoin(
        target_metric,
        source_metric,
        how="inner",
        predicate="intersects",
    )
    if joined.empty:
        return gpd.GeoDataFrame(columns=target_gdf.columns, crs=target_gdf.crs)

    target_areas = target_metric.geometry.area
    acc: dict[Any, float] = {}

    for target_idx, row in tqdm(
        joined.iterrows(),
        total=len(joined),
        desc="Area-weighted resampling",
        unit=" cells",
    ):
        try:
            source_idx = row["index_right"]
            target_geom = target_metric.at[target_idx, "geometry"]
            source_geom = source_metric.at[source_idx, "geometry"]
            target_area = target_areas.at[target_idx]
            if target_area == 0:
                continue
            intersection = _safe_intersection(target_geom, source_geom)
            if intersection.is_empty:
                continue
            proportion = intersection.area / target_area
            source_val = float(source_metric.at[source_idx, resample_col])
            acc[target_idx] = acc.get(target_idx, 0.0) + source_val * proportion
        except Exception:
            continue

    if not acc:
        return gpd.GeoDataFrame(columns=target_gdf.columns, crs=target_gdf.crs)

    out_rows: list[dict] = []
    for target_idx, resampled_value in acc.items():
        rec = target_gdf.loc[target_idx].to_dict()
        rec[resample_col] = round(resampled_value, 3)
        out_rows.append(rec)

    return gpd.GeoDataFrame(out_rows, crs=target_gdf.crs)


def _resampling_nearest(
    source_gdf: gpd.GeoDataFrame,
    target_gdf: gpd.GeoDataFrame,
    resample_col: str,
) -> gpd.GeoDataFrame:
    if source_gdf.empty or target_gdf.empty:
        return gpd.GeoDataFrame(columns=target_gdf.columns, crs=target_gdf.crs)

    for val in source_gdf[resample_col]:
        if not isinstance(val, Number):
            raise TypeError(
                f"Non-numeric value found in <{resample_col}>. "
                "Resampled field calculation failed."
            )

    source_metric, target_metric = _reproject_for_metric(source_gdf, target_gdf)

    joined = gpd.sjoin(
        target_metric,
        source_metric,
        how="inner",
        predicate="intersects",
    )
    if joined.empty:
        return gpd.GeoDataFrame(columns=target_gdf.columns, crs=target_gdf.crs)

    target_indices = joined.index.unique()

    source_pts = gpd.GeoDataFrame(
        {resample_col: source_metric[resample_col].values},
        geometry=source_metric.geometry.centroid,
        crs=source_metric.crs,
    )
    target_pts = gpd.GeoDataFrame(
        index=target_indices,
        geometry=target_metric.loc[target_indices].geometry.centroid,
        crs=target_metric.crs,
    )
    nearest_joined = gpd.sjoin_nearest(
        target_pts,
        source_pts,
        how="left",
    )

    out_rows: list[dict] = []
    for target_idx in tqdm(
        target_indices,
        total=len(target_indices),
        desc="Nearest neighbor resampling",
        unit=" cells",
    ):
        rec = target_gdf.loc[target_idx].to_dict()
        resampled_value = float(nearest_joined.loc[target_idx, resample_col])
        rec[resample_col] = round(resampled_value, 3)
        out_rows.append(rec)

    return gpd.GeoDataFrame(out_rows, crs=target_gdf.crs)


def resampling(
    source_gdf: gpd.GeoDataFrame,
    target_gdf: gpd.GeoDataFrame,
    resample_col: str,
    method: str = "nearest",
) -> gpd.GeoDataFrame:
    """
    Transfer ``resample_col`` from source cells onto target cells.

    Parameters
    ----------
    method
        ``\"area_weighted\"`` (default) — for each target geometry::

            value = sum(source_value * area(intersection) / area(target_cell))

        Only target rows with at least one intersecting source cell are returned.

        ``\"nearest\"`` — assign the value of the source cell whose centroid is
        closest to the target centroid (``geopandas.sjoin_nearest``).
        Only target rows that intersect at least one source cell are returned.
    """
    if resample_col not in source_gdf.columns:
        raise ValueError(
            f"There is no <{resample_col}> column in the source GeoDataFrame."
        )

    norm = method.strip().lower().replace("-", "_")
    if norm in ("area_weighted", "area"):
        return _resampling_area_weighted(source_gdf, target_gdf, resample_col)
    if norm in ("nearest", "nn", "nearest_neighbour", "nearest_neighbor"):
        return _resampling_nearest(source_gdf, target_gdf, resample_col)

    raise ValueError(
        f"Unsupported resampling method {method!r}; use 'area_weighted' or 'nearest'."
    )


def dggsresample(
    source_dggs: Union[gpd.GeoDataFrame, Any],
    dggs_from: str,
    dggs_to: str,
    resolution: Optional[int] = None,
    dggs_col: Optional[str] = None,
    resample_col: Optional[str] = None,
    output_format: str = "gpd",
    output_name: Optional[str] = None,
    method: str = "area_weighted",
    *,
    fix_antimeridian: Optional[str] = None,
    a5_options: Optional[dict] = None,
    split_antimeridian: bool = False,
    aggregate: bool = False,
    dggrid_options: Optional[dict] = None,
) -> Union[gpd.GeoDataFrame, str, dict, list, None]:
    """
    High-level resample: optional automatic resolution (-1), target grid
    generation, then optional attribute transfer.

    Parameters
    ----------
    source_dggs
        Source DGGS cells (GeoDataFrame, GeoJSON dict, feature list, or file path)
        with polygon geometry, DGGS id column, and optional numeric field.
    dggs_from,  dggs_to
        DGGS identifiers (e.g. ``\"h3\"``, ``\"s2\"``). DGGAL uses either the
        short type or the prefixed id column style: ``\"dggal_gnosis\"``,
        ``\"dggal_isea3h\"``, ``\"dggal_rtea7h\"``, etc. (see ``DGGAL_TYPES``).
        Bare ``\"rhealpix\"`` is vgrid's native rHEALPix; use ``\"dggal_rhealpix\"``
        for DGGAL's rHEALPix.
        DGGRID uses ``\"dggrid_ISEA4T\"``, ``\"dggrid_ISEA7H\"``, etc., or bare
        keys exactly as in ``DGGRID_TYPES`` (e.g. ``\"ISEA7H\"``). Native
        ``\"isea4t\"`` remains the Windows EAGGR grid, not DGGRID.
    resolution
        Target resolution (int), or ``None`` / ``-1`` to pick the nearest level
        by mean cell area.
    dggs_col
        Column holding source cell ids; defaults to ``dggs_from``. For DGGRID
        sources, if a canonical ``dggrid_<type>`` column exists (lowercase type,
        as from vgrid binning/conversion), it is used when the default name is
        missing.
    resample_col
        Numeric column on the source layer to redistribute; if omitted, only the
        target grid is returned.
    method
        ``\"area_weighted\"`` (default) or ``\"nearest\"`` (nearest source
        centroid via ``geopandas.sjoin_nearest``).
    output_format
        Output format; see :func:`~vgrid.utils.io.convert_to_output_format`.
    output_name
        Base name or filename for file-based outputs.
    aggregate
        When ``to_dggs`` is a DGGRID type (``dggrid_*``), dissolve split geometries
        by ``global_id`` after generation (requires antimeridian splitting).
    dggrid_options
        Optional dict passed to DGGRID ``grid_cell_polygons_for_extent`` (e.g.
        ``{\"densification\": 2}``). Only used when ``to_dggs`` is a DGGRID type.
    """
    if dggs_col is None:
        dggs_col = dggs_from

    source_gdf = process_input_data_resample(
        source_dggs,
        dggs_col=dggs_col,
        resample_col=resample_col,
    )

    dj = _dggrid_short_type(dggs_from)
    if dj is not None:
        dj = validate_dggrid_type(dj)
        canon = f"dggrid_{dj.lower()}"
        if canon in source_gdf.columns and dggs_col not in source_gdf.columns:
            dggs_col = canon
        elif canon in source_gdf.columns and dggs_col == dggs_from:
            dggs_col = canon

    if dggs_col not in source_gdf.columns:
        raise ValueError(f"Missing '{dggs_col}' in input data.")

    if resolution is None or resolution == -1:
        res = get_nearest_resolution(source_gdf, dggs_from, dggs_to, dggs_col)
    else:
        res = resolution

    target_gdf = generate_grid(
        source_gdf,
        dggs_to,
        int(res),
        fix_antimeridian=fix_antimeridian,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        dggrid_options=dggrid_options,
        a5_options=a5_options,  # for A5 grid generation
    )

    if resample_col:
        target_gdf = resampling(source_gdf, target_gdf, resample_col, method=method)

    if output_name is None and output_format in OUTPUT_FORMATS:
        if isinstance(source_dggs, str):
            base = os.path.splitext(os.path.basename(source_dggs))[0]
            output_name = f"{base}_{dggs_from}_to_{dggs_to}"
        else:
            output_name = f"{dggs_from}_to_{dggs_to}"
    return convert_to_output_format(target_gdf, output_format, output_name)


def dggsresample_cli():
    """Command-line interface for :func:`dggsresample`."""
    parser = argparse.ArgumentParser(
        description=(
            "Resample DGGS cells from one grid to another "
            "(area-weighted overlap or nearest-neighbour transfer)"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Source DGGS layer (vector file path, URL, GeoJSON, or GeoDataFrame path)",
    )
    parser.add_argument(
        "--from_dggs",
        "--dggs_from",
        dest="dggs_from",
        type=str,
        required=True,
        help=(
            "Source DGGS type (e.g. h3, s2, rhealpix, dggal_gnosis, dggrid_ISEA7H). "
            "Bare rhealpix is vgrid native; use dggal_rhealpix for DGGAL."
        ),
    )
    parser.add_argument(
        "--to_dggs",
        "--dggs_to",
        dest="dggs_to",
        type=str,
        required=True,
        help="Target DGGS type (same naming as --from_dggs)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=None,
        help="Target resolution; omit or pass -1 to match mean source cell area",
    )
    parser.add_argument(
        "-dggs_col",
        "--dggs_col",
        type=str,
        default=None,
        help="Source cell id column (default: same as --from_dggs)",
    )
    parser.add_argument(
        "-resample_col",
        "--resample_col",
        type=str,
        default=None,
        help="Numeric attribute to transfer; omit to return target grid only",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="area_weighted",
        choices=[
            "area_weighted",
            "area",
            "nearest",
            "nn",
            "nearest_neighbour",
            "nearest_neighbor",
        ],
        help="Resampling method (default: area_weighted)",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format (default: gpd)",
    )
    parser.add_argument(
        "-o",
        "--output_name",
        type=str,
        default=None,
        help="Output base name or file path for file-based formats",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method for supported grids",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Split geometries at the antimeridian",
    )
    parser.add_argument(
        "-aggregate",
        "--aggregate",
        action="store_true",
        default=False,
        help="For DGGRID targets: dissolve split geometries by global_id",
    )
    parser.add_argument(
        "--dggrid_options",
        type=str,
        default=None,
        help="JSON options for DGGRID grid generation (e.g. '{\"densification\": 2}')",
    )

    parser.add_argument(
        "--a5_options",
        "--a5_options",
        type=str,
        default=None,
        help="JSON options for A5 grid generation (e.g. '{\"segments\": 1000}')",
    )

    args = parser.parse_args()

    a5_options = None
    if args.a5_options:
        try:
            a5_options = json.loads(args.a5_options)
        except json.JSONDecodeError as exc:
            print(f"Error: Invalid JSON in a5_options: {exc}", file=sys.stderr)
            sys.exit(1)

    dggrid_options = None
    if args.dggrid_options:
        try:
            dggrid_options = json.loads(args.dggrid_options)
        except json.JSONDecodeError as exc:
            print(f"Error: Invalid JSON in dggrid_options: {exc}", file=sys.stderr)
            sys.exit(1)

    try:
        result = dggsresample(
            source_dggs=args.input,
            dggs_from=args.dggs_from,
            dggs_to=args.dggs_to,
            resolution=args.resolution,
            dggs_col=args.dggs_col,
            resample_col=args.resample_col,
            output_format=args.output_format,
            output_name=args.output_name,
            method=args.method,
            fix_antimeridian=args.fix_antimeridian,
            split_antimeridian=args.split_antimeridian,
            aggregate=args.aggregate,
            dggrid_options=dggrid_options,
            a5_options=a5_options,  # for A5 grid generation
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    dggsresample_cli()

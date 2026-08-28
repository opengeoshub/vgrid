"""Shared valuation helpers. Tolerances are the spec; tests assert against them."""

from __future__ import annotations

import importlib
import platform
from typing import Any

import numpy as np
import pytest
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry

from vgrid.utils.constants import (
    AUTHALIC_AREA,
    DGGAL_TYPES,
    DGGS_TYPES,
    GEOREF_RESOLUTION_DEGREES,
)
from vgrid.utils.io import gars_num_cells

N_SAMPLES = 1000
N_NATIVE = 200
RNG_SEED = 42
RES = 4
WORLD = [-180.0, -90.0, 180.0, 90.0]
WINDOWS_ONLY = frozenset()
ENUM_MAX_CELLS = 2000

NATIVE_DGGS = [
    name
    for name in DGGS_TYPES
    if name not in {"digipin", "ease", "isea3h", "isea4t", "mgrs"}
]
DGGAL_DGGS = list(DGGAL_TYPES.keys())

# ---------------------------------------------------------------------------
# Tolerances (stated; tests fail if exceeded)
# ---------------------------------------------------------------------------
# check                  | unit                         | value
# roundtrip contain      | fail rate over N_SAMPLES     | 0
# roundtrip buffer       | deg (planar fallback)        | 1e-6
# cell count             | absolute integer             | 0
# tessellated area       | relative |sum-A|/A           | 1e-9
# native cell id         | exact match                  | 0
# native geometry        | Hausdorff distance (deg)     | 1e-9

ROUNDTRIP_FAIL_RATE = 0.0
ROUNDTRIP_BUFFER_DEG = 1e-4
# Planar Shapely vs geodesic cell edges (antimeridian / poles).
ROUNDTRIP_FAIL_RATE_BY = {
    "h3": 5e-3,
    "s2": 5e-2,
    "a5": 5e-3,
    "healpix": 1e-2,
    "rhealpix": 2e-2,
    "qtm": 8e-2,
    "tilecode": 1e-2,
    "quadkey": 1e-2,
}
COUNT_ABS_TOL = 0
AREA_REL_TOL = 1e-9
AREA_REL_TOL_BY = {
    "dggal_rhealpix": 1e-3,  # polar dart cells
}

NATIVE_ID_ABS_TOL = 0
NATIVE_GEOM_HAUSDORFF_DEG = 1e-9

TOLERANCE_TABLE = [
    ("roundtrip_fail_rate", "fraction of N_SAMPLES (graticule)", ROUNDTRIP_FAIL_RATE),
    *[
        (f"roundtrip_fail_rate[{name}]", "fraction of N_SAMPLES", rate)
        for name, rate in ROUNDTRIP_FAIL_RATE_BY.items()
    ],
    ("roundtrip_buffer", "deg", ROUNDTRIP_BUFFER_DEG),
    ("cell_count", "absolute", COUNT_ABS_TOL),
    ("tessellated_area", "relative |sum-A|/AUTHALIC_AREA", AREA_REL_TOL),
    ("tessellated_area[dggal_rhealpix]", "relative vs AUTHALIC_AREA", AREA_REL_TOL_BY["dggal_rhealpix"]),
    ("native_cell_id", "absolute (exact)", NATIVE_ID_ABS_TOL),
    ("native_geometry", "Hausdorff deg", NATIVE_GEOM_HAUSDORFF_DEG),
]

TESSELLATION_RES = {
    "h3": 0,
    "s2": 0,
    "a5": 0,
    "healpix": 0,
    "rhealpix": 0,
    "qtm": 1,
    "olc": 2,
    "geohash": 1,
    "georef": 0,
    "tilecode": 1,
    "quadkey": 1,
    "maidenhead": 1,
}

TESSELLATION_SKIP = frozenset({"gars"})


def format_tolerance_table() -> str:
    headers = ("check", "unit", "tolerance")
    rows = [(name, unit, _fmt_tol(val)) for name, unit, val in TOLERANCE_TABLE]
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))

    def fmt(row) -> str:
        return " | ".join(str(c).ljust(widths[i]) for i, c in enumerate(row))

    sep = "-+-".join("-" * w for w in widths)
    return "\n".join([fmt(headers), sep, *[fmt(r) for r in rows]])


def _fmt_tol(val) -> str:
    if isinstance(val, str):
        return val
    if val == 0:
        return "0"
    return f"{val:.0e}"


def skip_windows_only(name: str) -> None:
    if name in WINDOWS_ONLY and platform.system() != "Windows":
        pytest.skip(f"{name} is only supported on Windows")


def sample_res(name: str) -> int:
    spec = DGGS_TYPES[name]
    lo, hi = spec["min_res"], spec["max_res"]
    return min(max(RES, lo), hi)


def dggal_sample_res(dggs_type: str) -> int:
    spec = DGGAL_TYPES[dggs_type]
    lo, hi = spec["min_res"], spec["max_res"]
    return min(max(RES, lo), hi)


def points_on_sphere(n: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform on the sphere. Returns (n, 2) as (lat, lon) degrees."""
    vec = rng.normal(size=(n, 3))
    vec /= np.linalg.norm(vec, axis=1, keepdims=True)
    lat = np.degrees(np.arcsin(np.clip(vec[:, 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(vec[:, 1], vec[:, 0]))
    return np.column_stack([lat, lon])


def sample_points(name: str, n: int, rng: np.random.Generator) -> np.ndarray:
    return points_on_sphere(n, rng)


def latlon_fn(name: str):
    mod = importlib.import_module("vgrid.conversion.latlon2dggs")
    return getattr(mod, f"latlon2{name}")


def geo_fn(name: str):
    mod = importlib.import_module(f"vgrid.conversion.dggs2geo.{name}2geo")
    return getattr(mod, f"{name}2geo")


def as_geometry(result: Any) -> BaseGeometry:
    if result is None:
        pytest.fail("2geo returned None")
    if isinstance(result, str):
        pytest.fail(f"2geo returned {result!r}")
    if isinstance(result, (list, tuple)):
        assert result, "2geo returned an empty list"
        result = result[0]
    if not isinstance(result, BaseGeometry):
        pytest.fail(f"2geo returned {type(result)!r}")
    return result


def point_in_cell(geom: BaseGeometry, lat: float, lon: float) -> bool:
    geom = as_geometry(geom)
    if geom.is_empty:
        return False
    candidates = [lon, lon - 360.0, lon + 360.0]
    buf = None
    for wrap in candidates:
        pt = Point(wrap, lat)
        if geom.covers(pt):
            return True
        if isinstance(geom, (Polygon, MultiPolygon)):
            if buf is None:
                buf = geom.buffer(ROUNDTRIP_BUFFER_DEG)
            if buf.covers(pt):
                return True
    return False


def analytic_count(name: str, res: int) -> int:
    if name == "h3":
        return 2 + 120 * (7**res)
    if name == "s2":
        return 6 * (4**res)
    if name == "a5":
        return 12 if res == 0 else 12 * 5 * (4 ** (res - 1))
    if name == "healpix":
        return 12 * (4**res)
    if name == "rhealpix":
        return 6 * (9**res)
    if name == "qtm":
        return 8 * (4 ** (res - 1))
    if name == "olc":
        if res <= 10:
            return 162 * (400 ** ((res // 2) - 1))
        return 162 * (400**4) * (20 ** (res - 10))
    if name == "geohash":
        return 32**res
    if name == "georef":
        deg = GEOREF_RESOLUTION_DEGREES[res]
        return int(round(360.0 / deg)) * int(round(180.0 / deg))
    if name == "tilecode":
        return 4**res
    if name == "quadkey":
        return 4**res
    if name == "maidenhead":
        from vgrid.dggs import maidenhead

        return maidenhead.num_cells(res)
    if name == "gars":
        return gars_num_cells(res)
    raise KeyError(name)


def native_count(name: str, res: int) -> int | None:
    """Library count when the backend exposes one; else None."""
    if name == "h3":
        import h3

        return int(h3.get_num_cells(res))
    if name == "a5":
        from a5.core.cell_info import get_num_cells

        return int(get_num_cells(res))
    if name == "healpix":
        from vgrid.dggs.healpix import nside2npix, order2nside

        return int(nside2npix(order2nside(res)))
    if name == "rhealpix":
        from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
        from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

        dggs = RHEALPixDGGS(ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=3, N_side=3)
        return int(dggs.num_cells(res))
    return None


def _dggal_dggrs(dggs_type: str):
    import dggal

    app = dggal.Application(appGlobals=dggal.__dict__)
    dggal.pydggal_setup(app)
    return getattr(dggal, DGGAL_TYPES[dggs_type]["class_name"])()


def dggal_native_count(dggs_type: str, res: int) -> int:
    return int(_dggal_dggrs(dggs_type).countZones(res))


def area_rel_tol(name: str) -> float:
    return AREA_REL_TOL_BY.get(name, AREA_REL_TOL)


def expected_area(name: str) -> float:
    return float(AUTHALIC_AREA)


def count_res_range(name: str):
    spec = DGGS_TYPES[name]
    lo, hi = spec["min_res"], spec["max_res"]
    if name == "olc":
        from vgrid.utils.io import olc_resolutions

        return [r for r in olc_resolutions if r <= lo + 8][:6]
    return range(lo, min(lo + 6, hi + 1))


def generate_grid(name: str, res: int):
    """Full-globe GeoDataFrame at `res`. verbose off."""
    skip_windows_only(name)
    if name == "h3":
        from vgrid.generator.h3grid import h3_grid

        return h3_grid(res, verbose=False)
    if name == "s2":
        from vgrid.generator.s2grid import s2_grid

        return s2_grid(res, WORLD, verbose=False)
    if name == "a5":
        from vgrid.generator.a5grid import a5_grid

        return a5_grid(res, WORLD, verbose=False)
    if name == "healpix":
        from vgrid.generator.healpixgrid import healpix_grid

        return healpix_grid(res, verbose=False)
    if name == "rhealpix":
        from vgrid.generator.rhealpixgrid import rhealpix_grid

        return rhealpix_grid(res, verbose=False)
    if name == "qtm":
        from vgrid.generator.qtmgrid import qtm_grid

        return qtm_grid(res, verbose=False)
    if name == "olc":
        from vgrid.generator.olcgrid import olc_grid

        return olc_grid(res, verbose=False)
    if name == "geohash":
        from vgrid.generator.geohashgrid import geohash_grid

        return geohash_grid(res, verbose=False)
    if name == "georef":
        from vgrid.generator.georefgrid import georef_grid

        return georef_grid(res, verbose=False)
    if name == "tilecode":
        from vgrid.generator.tilecodegrid import tilecode_grid

        return tilecode_grid(res, WORLD, verbose=False)
    if name == "quadkey":
        from vgrid.generator.quadkeygrid import quadkey_grid

        return quadkey_grid(res, WORLD, verbose=False)
    if name == "maidenhead":
        from vgrid.generator.maidenheadgrid import maidenhead_grid

        return maidenhead_grid(res, verbose=False)
    raise KeyError(name)


def generate_ids(name: str, res: int) -> list:
    skip_windows_only(name)
    if name == "h3":
        from vgrid.generator.h3grid import h3_grid_ids

        return list(h3_grid_ids(res, verbose=False))
    if name == "s2":
        from vgrid.generator.s2grid import s2_grid_ids

        return list(s2_grid_ids(res, WORLD))
    if name == "a5":
        from vgrid.generator.a5grid import a5_grid_ids

        return list(a5_grid_ids(res, WORLD))
    if name == "healpix":
        from vgrid.generator.healpixgrid import healpix_grid_ids

        return list(healpix_grid_ids(res))
    if name == "rhealpix":
        from vgrid.generator.rhealpixgrid import rhealpix_grid_ids

        return list(rhealpix_grid_ids(res))
    if name == "qtm":
        from vgrid.generator.qtmgrid import qtm_grid_ids

        return list(qtm_grid_ids(res, verbose=False))
    if name == "olc":
        from vgrid.generator.olcgrid import olc_grid_ids

        return list(olc_grid_ids(res, verbose=False))
    if name == "geohash":
        from vgrid.generator.geohashgrid import geohash_grid_ids

        return list(geohash_grid_ids(res))
    if name == "georef":
        from vgrid.generator.georefgrid import georef_grid_ids

        return list(georef_grid_ids(res, verbose=False))
    if name == "tilecode":
        from vgrid.generator.tilecodegrid import tilecode_grid_ids

        return list(tilecode_grid_ids(res))
    if name == "quadkey":
        from vgrid.generator.quadkeygrid import quadkey_grid_ids

        return list(quadkey_grid_ids(res))
    if name == "maidenhead":
        from vgrid.generator.maidenheadgrid import maidenhead_grid_ids

        return list(maidenhead_grid_ids(res, verbose=False))
    return None


def hausdorff_deg(a: BaseGeometry, b: BaseGeometry) -> float:
    return float(a.hausdorff_distance(b))

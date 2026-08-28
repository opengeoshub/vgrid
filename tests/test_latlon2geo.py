"""Round-trip tests for latlon2<dggs> and <dggs>2geo (Python API and CLI).

Sample point is from docs/notebooks/01_h3.ipynb.
"""

from __future__ import annotations

import importlib
import sys
from unittest.mock import patch

import pytest
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry

from vgrid.utils.constants import DGGAL_TYPES, DGGS_TYPES

LAT = 10.775276
LON = 106.706797
RES = 4

NATIVE_DGGS = [
    name
    for name in DGGS_TYPES
    if name not in {"digipin", "ease", "isea3h", "isea4t", "mgrs"}
]
DGGAL_DGGS = list(DGGAL_TYPES.keys())


def _import_latlon2dggs():
    try:
        return importlib.import_module("vgrid.conversion.latlon2dggs")
    except ImportError as exc:
        pytest.skip(f"Cannot import latlon2dggs: {exc}")


def _latlon_fn(name: str):
    return getattr(_import_latlon2dggs(), f"latlon2{name}")


def _latlon_cli(name: str):
    return getattr(_import_latlon2dggs(), f"latlon2{name}_cli")


def _geo_mod(name: str):
    try:
        return importlib.import_module(f"vgrid.conversion.dggs2geo.{name}2geo")
    except ImportError as exc:
        pytest.skip(f"Cannot import {name}2geo: {exc}")


def _geo_fn(name: str):
    return getattr(_geo_mod(name), f"{name}2geo")


def _geo_cli(name: str):
    return getattr(_geo_mod(name), f"{name}2geo_cli")


def _as_geometry(result) -> BaseGeometry:
    if result is None:
        pytest.fail("2geo returned None")
    if isinstance(result, (list, tuple)):
        assert result, "2geo returned an empty list"
        result = result[0]
    if not isinstance(result, BaseGeometry):
        pytest.fail(f"2geo returned {type(result)!r}, expected a Shapely geometry")
    return result


def _assert_covers_sample(geom, lat: float, lon: float) -> None:
    geom = _as_geometry(geom)
    assert not geom.is_empty
    point = Point(lon, lat)
    covers = geom.covers(point)
    if not covers and isinstance(geom, (Polygon, MultiPolygon)):
        covers = geom.buffer(1e-6).covers(point)
    assert covers, f"cell geometry does not cover ({lat}, {lon})"


def _run_cli(func, argv: list[str], capsys=None) -> str:
    with patch.object(sys, "argv", argv):
        func()
    if capsys is None:
        return ""
    return capsys.readouterr().out.strip()


def _cell_id_native(name: str):
    return _latlon_fn(name)(LAT, LON, RES)


def _require_in_coverage(name: str, cell_id):
    if str(cell_id) == "Out of Bound":
        pytest.skip(f"{name} does not cover sample point ({LAT}, {LON})")
    return cell_id


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_latlon2_python(dggs):
    cell_id = _cell_id_native(dggs)
    assert cell_id is not None
    assert str(cell_id) != ""


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_dggs2geo_python(dggs):
    cell_id = _require_in_coverage(dggs, _cell_id_native(dggs))
    geom = _geo_fn(dggs)(cell_id)
    _assert_covers_sample(geom, LAT, LON)


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_latlon2_cli(dggs, capsys):
    expected = str(_latlon_fn(dggs)(LAT, LON, RES))
    stdout = _run_cli(
        _latlon_cli(dggs),
        [f"latlon2{dggs}", str(LAT), str(LON), str(RES)],
        capsys,
    )
    assert stdout == expected


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_dggs2geo_cli(dggs):
    cell_id = _require_in_coverage(dggs, _cell_id_native(dggs))
    with patch.object(sys, "argv", [f"{dggs}2geo", str(cell_id)]):
        geom = _geo_cli(dggs)()
    _assert_covers_sample(geom, LAT, LON)


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_latlon2dggal_python(dggs_type):
    latlon = _import_latlon2dggs()
    zone_id = latlon.latlon2dggal(dggs_type, LAT, LON, RES)
    assert zone_id
    geom = _geo_fn("dggal")(dggs_type, zone_id)
    _assert_covers_sample(geom, LAT, LON)


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_latlon2dggal_cli(dggs_type, capsys):
    latlon = _import_latlon2dggs()
    expected = str(latlon.latlon2dggal(dggs_type, LAT, LON, RES))
    stdout = _run_cli(
        latlon.latlon2dggal_cli,
        ["latlon2dggal", dggs_type, str(LAT), str(LON), str(RES)],
        capsys,
    )
    assert stdout == expected


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_dggal2geo_cli(dggs_type):
    latlon = _import_latlon2dggs()
    zone_id = latlon.latlon2dggal(dggs_type, LAT, LON, RES)
    with patch.object(sys, "argv", ["dggal2geo", dggs_type, str(zone_id)]):
        geom = _geo_cli("dggal")()
    _assert_covers_sample(geom, LAT, LON)

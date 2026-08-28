"""Point → cell → geometry → point stays inside the cell.

Sample: N_SAMPLES points drawn uniformly on the sphere.
"""

from __future__ import annotations

import numpy as np
import pytest

from _common import (
    DGGAL_DGGS,
    N_SAMPLES,
    NATIVE_DGGS,
    RNG_SEED,
    ROUNDTRIP_FAIL_RATE,
    dggal_sample_res,
    geo_fn,
    latlon_fn,
    point_in_cell,
    sample_points,
    sample_res,
    skip_windows_only,
)


def _roundtrip(name, latlon, to_geo, pts, res, extra_args=()):
    misses = 0
    n = 0
    for lat, lon in pts:
        cell_id = latlon(*extra_args, lat, lon, res) if extra_args else latlon(lat, lon, res)
        if cell_id is None or str(cell_id) in {"", "None", "Out of Bound"}:
            continue
        n += 1
        geom = to_geo(cell_id) if not extra_args else to_geo(*extra_args, cell_id)
        if not point_in_cell(geom, float(lat), float(lon)):
            misses += 1
    if n == 0:
        pytest.fail(f"{name}: no in-coverage samples")
    rate = misses / n
    assert rate <= ROUNDTRIP_FAIL_RATE, (
        f"{name} res={res}: {misses}/{n} points outside cell (rate={rate:.4g} > {ROUNDTRIP_FAIL_RATE})"
    )


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_roundtrip_native(dggs):
    skip_windows_only(dggs)
    rng = np.random.default_rng(RNG_SEED)
    pts = sample_points(dggs, N_SAMPLES, rng)
    res = sample_res(dggs)
    _roundtrip(dggs, latlon_fn(dggs), geo_fn(dggs), pts, res)


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_roundtrip_dggal(dggs_type):
    from vgrid.conversion.latlon2dggs import latlon2dggal
    from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo

    rng = np.random.default_rng(RNG_SEED)
    pts = sample_points("dggal", N_SAMPLES, rng)
    res = dggal_sample_res(dggs_type)

    def latlon(lat, lon, resolution):
        return latlon2dggal(dggs_type, lat, lon, resolution)

    def to_geo(cell_id):
        return dggal2geo(dggs_type, cell_id)

    _roundtrip(f"dggal_{dggs_type}", latlon, to_geo, pts, res)

"""Sum of tessellated cell areas vs authalic ellipsoid area."""

from __future__ import annotations

import pytest

from _common import (
    DGGAL_DGGS,
    TESSELLATION_RES,
    TESSELLATION_SKIP,
    area_rel_tol,
    expected_area,
    generate_grid,
    skip_windows_only,
)
from vgrid.utils.constants import AUTHALIC_AREA


def _rel_err(sum_area: float, expected: float) -> float:
    return abs(float(sum_area) - expected) / expected


@pytest.mark.parametrize("dggs", sorted(TESSELLATION_RES))
def test_tessellated_area(dggs):
    skip_windows_only(dggs)
    if dggs in TESSELLATION_SKIP:
        pytest.skip(f"{dggs} global tessellation is too large for this check")
    res = TESSELLATION_RES[dggs]
    gdf = generate_grid(dggs, res)
    err = _rel_err(gdf.cell_area.sum(), expected_area(dggs))
    tol = area_rel_tol(dggs)
    assert err <= tol, (
        f"{dggs} res={res}: |sum-A|/A = {err:.3e} > {tol:.0e}"
    )


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_dggal_tessellated_area(dggs_type):
    from vgrid.generator.dggalgen import dggalgen

    gdf = dggalgen(dggs_type, 0, verbose=False, output_format="gdf")
    err = _rel_err(gdf.cell_area.sum(), AUTHALIC_AREA)
    tol = area_rel_tol(f"dggal_{dggs_type}")
    assert err <= tol, (
        f"dggal_{dggs_type} res=0: |sum-A|/A = {err:.3e} > {tol:.0e}"
    )

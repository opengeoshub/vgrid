"""Cell counts per resolution match the analytic formula (and the native lib when it exposes a count)."""

from __future__ import annotations

import pytest

from _common import (
    COUNT_ABS_TOL,
    DGGAL_DGGS,
    ENUM_MAX_CELLS,
    NATIVE_DGGS,
    analytic_count,
    count_res_range,
    dggal_native_count,
    generate_ids,
    native_count,
    skip_windows_only,
)


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_count_matches_analytic(dggs):
    skip_windows_only(dggs)
    lib_n = None
    for res in count_res_range(dggs):
        expected = analytic_count(dggs, res)
        lib_n = native_count(dggs, res)
        if lib_n is not None:
            assert abs(lib_n - expected) <= COUNT_ABS_TOL, (
                f"{dggs} res={res}: native {lib_n} != analytic {expected}"
            )
        if expected <= ENUM_MAX_CELLS:
            ids = generate_ids(dggs, res)
            if ids is None:
                continue
            assert abs(len(ids) - expected) <= COUNT_ABS_TOL, (
                f"{dggs} res={res}: generated {len(ids)} != analytic {expected}"
            )


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_dggal_count_matches_native(dggs_type):
    from vgrid.generator.dggalgen import dggalgen

    for res in range(0, 3):
        expected = dggal_native_count(dggs_type, res)
        if expected > ENUM_MAX_CELLS:
            continue
        gdf = dggalgen(dggs_type, res, verbose=False, output_format="gdf")
        n = len(gdf)
        assert abs(n - expected) <= COUNT_ABS_TOL, (
            f"dggal_{dggs_type} res={res}: generated {n} != countZones {expected}"
        )

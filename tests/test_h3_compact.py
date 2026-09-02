"""h3_compact(depth=-1) should match native h3.compact_cells."""

from __future__ import annotations

import h3
import pandas as pd
import pytest

from vgrid.conversion.dggscompact.h3compact import h3_compact, h3compact

LAT = 10.775276
LON = 106.706797


def _assert_same_as_native(h3_ids):
    native = set(h3.compact_cells(h3_ids))
    ours = set(h3_compact(h3_ids, depth=-1, verbose=False))
    assert ours == native


def test_depth_minus_one_matches_native_full_children():
    parent = h3.latlng_to_cell(LAT, LON, 3)
    children = list(h3.cell_to_children(parent, 5))
    _assert_same_as_native(children)
    assert h3_compact(children, depth=-1, verbose=False) == [parent]


def test_depth_minus_one_matches_native_incomplete_children():
    parent = h3.latlng_to_cell(LAT, LON, 3)
    children = list(h3.cell_to_children(parent, 5))[:-1]
    _assert_same_as_native(children)


def test_depth_minus_one_matches_native_grid_disk():
    center = h3.latlng_to_cell(LAT, LON, 8)
    disk = list(h3.grid_disk(center, 4))
    _assert_same_as_native(disk)


@pytest.mark.parametrize("res", [1, 4, 8])
def test_depth_minus_one_matches_native_across_resolutions(res):
    center = h3.latlng_to_cell(LAT, LON, res)
    disk = list(h3.grid_disk(center, 3))
    _assert_same_as_native(disk)


def test_depth_minus_one_matches_native_pentagon_children():
    pentagon = next(c for c in h3.get_res0_cells() if h3.is_pentagon(c))
    children = list(h3.cell_to_children(pentagon, 2))
    assert h3.is_pentagon(pentagon)
    _assert_same_as_native(children)


def test_depth_minus_one_matches_native_already_compact():
    cell = h3.latlng_to_cell(LAT, LON, 6)
    neighbors = list(h3.grid_disk(cell, 1))
    _assert_same_as_native(neighbors)


def test_compact_agg_count_one_level():
    parent = h3.latlng_to_cell(LAT, LON, 4)
    children = list(h3.cell_to_children(parent))
    gdf = pd.DataFrame({"h3": children})
    out = h3compact(gdf, depth=1, agg="count", verbose=False)
    assert len(out) == 1
    assert out.iloc[0]["h3"] == parent
    assert int(out.iloc[0]["count"]) == len(children)


def test_compact_agg_count_two_levels():
    parent = h3.latlng_to_cell(LAT, LON, 3)
    children = list(h3.cell_to_children(parent, 5))
    gdf = pd.DataFrame({"h3": children})
    out = h3compact(gdf, depth=-1, agg="count", verbose=False)
    assert len(out) == 1
    assert out.iloc[0]["h3"] == parent
    assert int(out.iloc[0]["count"]) == len(children)


def test_compact_agg_mean():
    parent = h3.latlng_to_cell(LAT, LON, 4)
    children = list(h3.cell_to_children(parent))
    values = list(range(1, len(children) + 1))
    gdf = pd.DataFrame({"h3": children, "value": values})
    out = h3compact(gdf, depth=1, agg="mean", numeric_col="value", verbose=False)
    assert len(out) == 1
    assert out.iloc[0]["h3"] == parent
    assert out.iloc[0]["value_mean"] == pytest.approx(sum(values) / len(values))


def test_compact_agg_sum_incomplete_children():
    parent = h3.latlng_to_cell(LAT, LON, 4)
    children = list(h3.cell_to_children(parent))[:-1]
    values = [10] * len(children)
    gdf = pd.DataFrame({"h3": children, "value": values})
    out = h3compact(gdf, depth=1, agg="sum", numeric_col="value", verbose=False)
    assert len(out) == len(children)
    assert set(out["h3"]) == set(children)
    assert out["value_sum"].tolist() == values


def test_compact_agg_requires_numeric_col():
    parent = h3.latlng_to_cell(LAT, LON, 4)
    children = list(h3.cell_to_children(parent))
    gdf = pd.DataFrame({"h3": children})
    with pytest.raises(ValueError, match="numeric_col"):
        h3compact(gdf, depth=1, agg="mean", verbose=False)

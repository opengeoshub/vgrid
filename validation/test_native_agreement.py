"""vgrid wrappers vs the native library called directly. Cell IDs exact; geometry Hausdorff."""

from __future__ import annotations

import re

import numpy as np
import pytest
from shapely.geometry import Polygon

from _common import (
    DGGAL_DGGS,
    N_NATIVE,
    NATIVE_DGGS,
    NATIVE_GEOM_HAUSDORFF_DEG,
    NATIVE_ID_ABS_TOL,
    RNG_SEED,
    as_geometry,
    dggal_sample_res,
    geo_fn,
    hausdorff_deg,
    latlon_fn,
    sample_points,
    sample_res,
    skip_windows_only,
)
from vgrid.utils.constants import DGGAL_TYPES


def _native_cell(name: str, lat: float, lon: float, res: int):
    if name == "h3":
        import h3

        return h3.latlng_to_cell(lat, lon, res)
    if name == "s2":
        from vgrid.dggs import s2

        cell = s2.CellId.from_lat_lng(s2.LatLng.from_degrees(lat, lon)).parent(res)
        return s2.CellId.to_token(cell)
    if name == "a5":
        import a5

        return a5.u64_to_hex(a5.lonlat_to_cell((lon, lat), res))
    if name == "healpix":
        from vgrid.dggs.healpix import lonLat2PixNest, order2nside, orderpix2uniq

        return orderpix2uniq(res, lonLat2PixNest(order2nside(res), lon, lat))
    if name == "rhealpix":
        from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
        from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

        dggs = RHEALPixDGGS(
            ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=3, N_side=3
        )
        return str(dggs.cell_from_point(res, (lon, lat), plane=False))
    if name == "qtm":
        from vgrid.dggs import qtm

        return qtm.latlon_to_qtm_id(lat, lon, res)
    if name == "olc":
        from vgrid.dggs import olc

        return olc.encode(lat, lon, res)
    if name == "geohash":
        from vgrid.dggs import geohash

        return geohash.encode(lat, lon, res)
    if name == "georef":
        from vgrid.dggs import georef

        return georef.encode(lat, lon, res)
    if name == "tilecode":
        from vgrid.dggs import tilecode

        return tilecode.latlon2tilecode(lat, lon, res)
    if name == "quadkey":
        from vgrid.dggs import tilecode

        return tilecode.latlon2quadkey(lat, lon, res)
    if name == "maidenhead":
        from vgrid.dggs import maidenhead

        return maidenhead.toMaiden(lat, lon, res)
    if name == "gars":
        from gars_field.garsgrid import GARSGrid

        minutes = {1: 30, 2: 15, 3: 5, 4: 1}[res]
        return str(GARSGrid.from_latlon(lat, lon, minutes))
    raise KeyError(name)


def _native_poly(name: str, cell_id):
    if name == "h3":
        import h3

        return Polygon([(lng, lat) for lat, lng in h3.cell_to_boundary(cell_id)])
    if name == "s2":
        from vgrid.dggs import s2

        cell = s2.Cell(s2.CellId.from_token(cell_id))
        verts = []
        for i in range(4):
            ll = s2.LatLng.from_point(cell.get_vertex(i))
            verts.append((ll.lng().degrees, ll.lat().degrees))
        verts.append(verts[0])
        return Polygon(verts)
    if name == "a5":
        import a5

        return Polygon(a5.cell_to_boundary(a5.hex_to_u64(cell_id)))
    if name == "healpix":
        from vgrid.dggs.healpix import cornersNestLonLat, order2nside, uniq2orderpix

        order, ipix = uniq2orderpix(int(cell_id))
        return Polygon(cornersNestLonLat(order2nside(order), ipix))
    if name == "rhealpix":
        from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
        from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
        from vgrid.utils.geometry import rhealpix_cell_to_polygon

        dggs = RHEALPixDGGS(
            ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=3, N_side=3
        )
        uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
        return rhealpix_cell_to_polygon(dggs.cell(uids))
    if name == "qtm":
        from vgrid.dggs.qtm import constructGeometry, qtm_id_to_facet

        return constructGeometry(qtm_id_to_facet(cell_id))
    if name == "olc":
        from vgrid.dggs import olc

        coord = olc.decode(cell_id)
        return Polygon(
            [
                [coord.longitudeLo, coord.latitudeLo],
                [coord.longitudeHi, coord.latitudeLo],
                [coord.longitudeHi, coord.latitudeHi],
                [coord.longitudeLo, coord.latitudeHi],
                [coord.longitudeLo, coord.latitudeLo],
            ]
        )
    if name == "geohash":
        from vgrid.dggs import geohash

        bbox = geohash.bbox(cell_id)
        return Polygon(
            [
                [bbox["w"], bbox["s"]],
                [bbox["e"], bbox["s"]],
                [bbox["e"], bbox["n"]],
                [bbox["w"], bbox["n"]],
                [bbox["w"], bbox["s"]],
            ]
        )
    if name == "georef":
        from vgrid.dggs import georef

        _, _, min_lat, min_lon, max_lat, max_lon, _ = georef.georefcell(cell_id)
        return Polygon(
            [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat],
            ]
        )
    if name == "tilecode":
        from vgrid.dggs import mercantile

        match = re.match(r"z(\d+)x(\d+)y(\d+)", str(cell_id))
        z, x, y = map(int, match.groups())
        b = mercantile.bounds(x, y, z)
        return Polygon(
            [
                [b.west, b.south],
                [b.east, b.south],
                [b.east, b.north],
                [b.west, b.north],
                [b.west, b.south],
            ]
        )
    if name == "quadkey":
        from vgrid.dggs import mercantile

        tile = mercantile.quadkey_to_tile(cell_id)
        b = mercantile.bounds(tile.x, tile.y, tile.z)
        return Polygon(
            [
                [b.west, b.south],
                [b.east, b.south],
                [b.east, b.north],
                [b.west, b.north],
                [b.west, b.south],
            ]
        )
    if name == "maidenhead":
        from vgrid.dggs import maidenhead

        _, _, min_lat, min_lon, max_lat, max_lon = maidenhead.maidenGrid(cell_id)
        return Polygon(
            [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat],
            ]
        )
    if name == "gars":
        from gars_field import garsgrid

        poly = garsgrid.GARSGrid(cell_id).polygon
        return Polygon(list(poly.exterior.coords))
    return None


@pytest.mark.parametrize("dggs", NATIVE_DGGS)
def test_native_agreement(dggs):
    skip_windows_only(dggs)
    rng = np.random.default_rng(RNG_SEED)
    pts = sample_points(dggs, N_NATIVE, rng)
    res = sample_res(dggs)
    wrap = latlon_fn(dggs)
    to_geo = geo_fn(dggs)
    n_id = 0
    n_geom = 0
    for lat, lon in pts:
        wrapped = wrap(lat, lon, res)
        native = _native_cell(dggs, float(lat), float(lon), res)
        if wrapped is None or native is None or str(wrapped) == "Out of Bound":
            continue
        assert str(wrapped) == str(native), f"{dggs}: wrapper {wrapped!r} != native {native!r}"
        native_geom = _native_poly(dggs, wrapped)
        if native_geom is None:
            n_id += 1
            continue
        wrap_geom = as_geometry(to_geo(wrapped))
        dist = hausdorff_deg(wrap_geom, native_geom)
        n_id += 1
        n_geom += 1
        assert dist <= NATIVE_GEOM_HAUSDORFF_DEG, (
            f"{dggs} {wrapped}: Hausdorff {dist:.3e} deg > {NATIVE_GEOM_HAUSDORFF_DEG:.0e}"
        )
    assert n_id > 0, f"{dggs}: no in-coverage samples"
    assert n_geom > 0, f"{dggs}: no geometry comparisons"


@pytest.mark.parametrize("dggs_type", DGGAL_DGGS)
def test_dggal_native_agreement(dggs_type):
    import dggal
    from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo
    from vgrid.conversion.latlon2dggs import latlon2dggal

    app = dggal.Application(appGlobals=dggal.__dict__)
    dggal.pydggal_setup(app)
    dggrs = getattr(dggal, DGGAL_TYPES[dggs_type]["class_name"])()

    rng = np.random.default_rng(RNG_SEED)
    pts = sample_points("dggal", N_NATIVE, rng)
    res = dggal_sample_res(dggs_type)
    compared_geom = 0
    for lat, lon in pts:
        wrapped = latlon2dggal(dggs_type, float(lat), float(lon), res)
        zone = dggrs.getZoneFromWGS84Centroid(res, dggal.GeoPoint(float(lat), float(lon)))
        native = dggrs.getZoneTextID(zone)
        assert str(wrapped) == str(native), (
            f"dggal_{dggs_type}: wrapper {wrapped!r} != native {native!r}"
        )
        wrap_geom = as_geometry(dggal2geo(dggs_type, wrapped))
        z = dggrs.getZoneFromTextID(wrapped)
        vertices = dggrs.getZoneRefinedWGS84Vertices(z, 0)
        coords = [
            (float(vertices[i].lon), float(vertices[i].lat))
            for i in range(vertices.count)
        ]
        coords.append(coords[0])
        native_geom = Polygon(coords)
        dist = hausdorff_deg(wrap_geom, native_geom)
        compared_geom += 1
        assert dist <= NATIVE_GEOM_HAUSDORFF_DEG, (
            f"dggal_{dggs_type} {wrapped}: Hausdorff {dist:.3e} deg > {NATIVE_GEOM_HAUSDORFF_DEG:.0e}"
        )
    assert compared_geom > 0

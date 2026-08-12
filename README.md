<!-- PROJECT LOGO -->
<p align="center">
  <strong >Vgrid </strong> <br>
    <b><i> A unified framework for working with DGGS</i><b>
</p>
<p align="center">
  <img src="https://raw.githubusercontent.com/opengeoshub/vgridtools/main/images/readme/dggs.png">
</p>


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Vgrid User Guide</summary>  
  <ol>
   <li>
      <a href="#vgrid-introduction">Vgrid Introduction</a>     
    </li>
    <li>
      <a href="#vgrid-installation">Vgrid Installation</a>     
    </li>
    <li>
      <a href="#dggs-conversion">DGGS Conversion</a>
      <ul>
        <li><a href="#lat-lon-to-dggs">Lat lon to DGGS</a></li>
        <li><a href="#dggs-to-geojson">DGGS to GeoJSON</a></li>
        <li><a href="#vector-to-dggs">Vector to DGGS</a></li>
        <li><a href="#dggs-compact">DGGS Compact</a></li>
        <li><a href="#dggs-expand">DGGS Expand</a></li>
        <li><a href="#dggs-resample">DGGS Resample</a></li>
        <li><a href="#raster-to-dggs">Raster to DGGS</a></li>
        <li><a href="#dggs-binning">DGGS Binning</a></li>
      </ul>
    </li>
    <li><a href="#dggs-generator">DGGS Generator</a></li>
    <li><a href="#dggs-inspect">DGGS Inspect</a></li>
    <li><a href="#dggs-stats">DGGS Stats</a></li>
    <li><a href="#polyhedra-generator">Polyhedra Generator</a></li>
  </ol>
</details>

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opengeoshub/vgrid/blob/main)
[![image](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/opengeoshub/vgrid/HEAD?filepath=docs/notebooks)
[![image](https://studiolab.sagemaker.aws/studiolab.svg)](https://studiolab.sagemaker.aws/import/github/opengeoshub/vgrid/blob/main/docs/notebooks/00_intro.ipynb)
[![image](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://jupyterlite.gishub.vn/lab/index.html?path=notebooks/vgrid/00_intro.ipynb)[![PyPI version](https://badge.fury.io/py/vgrid.svg)](https://badge.fury.io/py/vgrid)
[![image](https://static.pepy.tech/badge/vgrid)](https://pepy.tech/project/vgrid)
[![DOI](https://zenodo.org/badge/1251638045.svg)](https://doi.org/10.5281/zenodo.21717224)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Vgrid Introduction

Vgrid DGGS is a python package for working with popular geodesic DGGS and graticule-based DGGS, including conversion, compact/expand, resampling between grids, binning, raster ingestion, and grid generation.

Vgrid supports a wide range of popular DGGSs, including H3, S2, A5, rHEALPix, DGGAL, DGGRID, Open-EAGGR ISEA4T, ISEA3H, EASE-DGGS, QTM, OLC, Geohash, GEOREF, MGRS, TileCode, Quadkey, Maidenhead, and GARS.

Vgrid supports popular input and output GIS formats, including CSV, GeoJSON, Pandas/GeoPandas, Shapefile, GeoPackage, and GeoParquet.

***The following Vgrid user guide is intended for use in a Command Line Interface (CLI) environment. Full Vgrid DGGS documentation is available at [vgrid document](https://vgrid.gishub.vn).***

To work with Vgrid DGGS directly in GeoPandas and Pandas, please use [vgridpandas](https://pypi.org/project/vgridpandas/). Full Vgridpandas DGGS documentation is available at [vgridpandas document](https://vgridpandas.gishub.vn)

To work with Vgrid DGGS in QGIS, install the [Vgrid Plugin](https://plugins.qgis.org/plugins/vgridtools/).

To visualize DGGS in Maplibre GL JS, try the [vgrid-maplibre](https://www.npmjs.com/package/vgrid-maplibre) library.

For an interactive demo, visit the [Vgrid Homepage](https://vgridhome.gishub.vn).

## Citation

If you use Vgrid DGGS in your work, please cite it. Vgrid is archived on [Zenodo](https://zenodo.org/)

[![DOI](https://zenodo.org/badge/1251638045.svg)](https://doi.org/10.5281/zenodo.21717224)

### Acknowledgements
Vgrid is built upon free and open-source software and would like to acknowledge the maintainers and contributors of the following projects, together with the many transitive dependencies that make them possible.

- [h3-py](https://github.com/uber/h3-py) by [Uber](https://github.com/uber).
- [s2sphere](https://github.com/sidewalklabs/s2sphere) by [Sidewalk Labs](https://github.com/sidewalklabs).
- [a5-py](https://github.com/felixpalmer/a5-py) by [Felix Palmer](https://github.com/felixpalmer) and [Thang Quach](https://github.com/thangqd).
- [rhealpixdggs-py](https://github.com/manaakiwhenua/rhealpixdggs-py) by [Manaaki Whenua – Landcare Research](https://github.com/manaakiwhenua).
- [open-eaggr](https://github.com/riskaware-ltd/open-eaggr) by [Riskaware](https://github.com/riskaware-ltd).
- [EASE-DGGS](https://github.com/GEMS-UMN/EASE-DGGS) by [GEMS-UMN](https://github.com/GEMS-UMN).
- [pydggal](https://github.com/ecere/pydggal) by [Jerome St-Louis](https://github.com/jerstlouis) from [Ecere](https://github.com/ecere).
- [DGGRID](https://github.com/sahrk/DGGRID) by [Kevin Sahr](https://github.com/sahrk).
- [dggrid4py](https://github.com/allixender/dggrid4py) by [Alex Kmoch](https://github.com/allixender).
- [QTM](https://github.com/opengeoshub/vgrid/blob/main/vgrid/dggs/qtm.py) by [Thang Quach](https://github.com/thangqd), with reference to [QTM](https://github.com/paulojraposo/QTM) by [Paulo Raposo](https://github.com/paulojraposo).
- [Lat Lon Tools QGIS Plugin](https://github.com/hamiltoncj/qgis-latlontools-plugin) by [Calvin Hamilton](https://github.com/hamiltoncj).
- [geohash](https://github.com/hkwi/python-geohash/blob/master/geohash.py) by [Hiroaki Kawai](https://github.com/hkwi).
- [gars-field](https://github.com/corteva/gars-field) by [Corteva Agriscience](https://github.com/corteva).
- [Tilecode & Quadkey](https://github.com/opengeoshub/vgrid/blob/main/vgrid/dggs/tilecode.py) by [Thang Quach](https://github.com/thangqd), utilizing [mercantile](https://github.com/mapbox/mercantile) by [Mapbox](https://github.com/mapbox).
- [antimeridian](https://www.gadom.ski/antimeridian) by [gadomski](https://github.com/gadomski/antimeridian).
- The DGGS Inspect feature in Vgrid is inspired by [Area and shape distortions in open-source discrete global grid systems](https://www.tandfonline.com/doi/full/10.1080/20964471.2022.2094926#abstract) by [Alex Kmoch](https://github.com/allixender) et al. (2022) ([resources](https://github.com/LandscapeGeoinformatics/dggs_comparisons)).
- The [Vgrid Document](https://vgrid.gishub.vn/) is inspired by [leafmap](https://leafmap.org/) developed by [Qiusheng Wu](https://github.com/giswqs) from [Open Geospatial Solutions](https://github.com/opengeos).

## Vgrid Installation
- Using pip:   
    ``` bash 
    pip install vgrid --upgrade
    ```
## DGGS Conversion

### Lat lon to DGGS

Convert lat, long in WGS84 CRS to DGGS cellID (H3, S2, A5, rHEALPix, DGGAL, DGGRID, Open-EAGGR ISEA4T and ISEA3H, EASE-DGGS, QTM, OLC, Geohash, GEOREF, MGRS, Tilecode, Quadkey, Maidenhead, GARS, DIGIPIN).

``` bash
> latlon2h3 10.775276 106.706797 13 # latlon2h3 <lat> <lon> <resolution> [0..15] 
> latlon2s2 10.775276 106.706797 21 # latlon2s2 <lat> <lon> <resolution> [0..30]
> latlon2a5 10.775276 106.706797 18 # latlon2a5 <lat> <lon> <resolution> [0..29]
> latlon2rhealpix 10.775276 106.706797 14 # latlon2rhealpix <lat> <lon> <resolution> [0..15]
> latlon2dggal 10.775276 106.706797 gnosis 5 # latlon2dggal <lat> <lon> <dggal_type> <resolution>
> latlon2dggrid 10.775276 106.706797 FULLER4D 10 # Linux only: latlon2dggrid <lat> <lon> <dggs_type> <resolution>
> latlon2isea4t 10.775276 106.706797 21 # latlon2isea4t <lat> <lon> <resolution> [0..39]
> latlon2isea3h 10.775276 106.706797 27 # latlon2isea3h <lat> <lon> <resolution> [0..40]
> latlon2ease 10.775276 106.706797 6 # latlon2ease <lat> <lon> <resolution> [0..6]
> latlon2qtm 10.775276 106.706797 11  # latlon2qtm <lat> <lon> <resolution> [1..24]
> latlon2olc 10.775276 106.706797 11 # latlon2olc <lat> <lon> <resolution> [2,4,6,8,10..15]
> latlon2geohash 10.775276 106.706797 9 # latlon2geohash <lat> <lon> <resolution> [1..10]
> latlon2georef 10.775276 106.706797 4 # latlon2georef <lat> <lon> <resolution> [0..10]
> latlon2mgrs 10.775276 106.706797 4 # latlon2mgrs <lat> <lon> <resolution> [0..5]
> latlon2tilecode 10.775276 106.706797 23 # latlon2tilecode <lat> <lon> <resolution> [0..29]
> latlon2quadkey 10.775276 106.706797 23 # latlon2quadkey <lat> <lon> <resolution> [0..29]
> latlon2maidenhead 10.775276 106.706797 4 # latlon2maidenhead <lat> <lon> <resolution> [1..4]
> latlon2gars 10.775276 106.706797 1 # latlon2gars <lat> <lon> <resolution> [1..4] (30,15,5,1 minutes)
> latlon2digipin 10.775276 106.706797 10 # latlon2digipin <lat> <lon> <resolution>
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/latlon2dggs.png">
</div>

### DGGS to GeoJSON
Convert DGGS cell IDs to GeoJSON.

``` bash
> h32geojson 8d65b56628e46bf 
> s22geojson 31752f45cc94 
> a52geojson 8e65b56628e0d07
> rhealpix2geojson R31260335553825
> dggal2geojson isea3h A4-0-A # dggal2geojson <dggal_type> <zone_id> [-split]
> dggrid2geojson 8478420 FULLER4D 10 # Linux only: dggrid2geojson <DGGRID code> <DGGS Type> <resolution>
# <DGGS Type>: SUPERFUND, PLANETRISK, ISEA3H, ISEA4H, ISEA4T, ISEA4D, ISEA43H, ISEA7H, IGEO7, FULLER3H, FULLER4H, FULLER4T, FULLER4D, FULLER43H, FULLER7H
> isea4t2geojson 13102313331320133331133 
> isea3h2geojson 1327916769,-55086 
> ease2geojson L6.165767.02.02.22.45.63.05
> qtm2geojson 42012323131
> olc2geojson 7P28QPG4+4P7
> geohash2geojson w3gvk1td8
> georef2geojson VGBL42404651
> mgrs2geojson 34TGK56063228
> gzd # Create Grid Zone Designators - used by MGRS
> tilecode2geojson z23x6680749y3941729
> quadkey2geojson 13223011131020212310000
> maidenhead2geojson OK30is46 
> gars2geojson 574JK1918
> digipin2geojson F3K
```

```geojson
{"type": "FeatureCollection", "features": [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [[[106.70789209957029, 10.77601109480422], [106.70745212203926, 10.777918172217145], [106.7055792173179, 10.778494886227104], [106.70414629547702, 10.777164500398262], [106.70458630070586, 10.775257408234628], [106.70645920007833, 10.774680716650128], [106.70789209957029, 10.77601109480422]]]}, "properties": {"h3": "8965b56628fffff", "resolution": 9, "center_lat": 10.7765878, "center_lon": 106.7060192, "avg_edge_len": 215.295, "cell_area": 120421.396}}]}
```

### Vector to DGGS
Convert vector layers (Point/MultiPoint, LineString/MultiLineString, Polygon/MultiPolygon) to DGGS cells.

Common `vector2*` options: `-i` / `--input`, `-r` / `--resolution`, `-p` / `--predicate` (`intersect`, `within`, `centroid_within`, `largest_overlap` for polygons), `-c` / `--compact`, `-t` / `--topology`, `-f` / `--output_format`. Antimeridian: `-fix` for H3, S2, rHEALPix, ISEA4T/ISEA3H; `-split` for A5, DGGAL, DGGRID; A5/DGGRID also support `-options` / `--options` (JSON).

``` bash
> vector2h3 -i polygon.geojson -r 11 -p intersect -c -fix split -f geojson # vector2h3 -i <input> -r <resolution[0..15]> [-p predicate] [-c] [-t] [-f output_format] [-fix method]
> vector2s2 -i polygon.geojson -r 18 -p intersect -c -fix split # vector2s2 -r <resolution[0..30]>
> vector2a5 -i polyline.geojson -r 18 -p intersect -split # vector2a5 -r <resolution[0..29]>; polylines use great-circle path cells
> vector2rhealpix -i polygon.geojson -r 11 -p largest_overlap -fix split # vector2rhealpix -r <resolution[0..15]>
> vector2dggal -i polygon.geojson -dggs gnosis -r 8 -p intersect -c # vector2dggal -i <input> -dggs <dggal_type> -r <resolution> [-p predicate] [-c compact] [-t topology]
> vector2dggrid -i polyline.geojson -t ISEA4H -r 17 -split -aggregate # Linux: vector2dggrid -i <input> -t <DGGS Type> -r <resolution> [-a address_type] [-split] [-aggregate] [-options JSON]
> vector2isea4t -i polygon.geojson -r 17 -p intersect # vector2isea4t -r <resolution[0..25]> (Windows)
> vector2isea3h -i polygon.geojson -r 18 -p intersect # vector2isea3h -r <resolution[0..40]> (Windows)
> vector2ease -i polygon.geojson -r 4 -p intersect # vector2ease -r <resolution[0..6]>
> vector2qtm -i polygon.geojson -r 18 -p intersect -c # vector2qtm -r <resolution[1..24]>
> vector2olc -i polygon.geojson -r 10 -p intersect # vector2olc -r <resolution[2,4,6,8,10..15]>
> vector2geohash -i polygon.geojson -r 7 -p intersect # vector2geohash -r <resolution[1..10]>
> vector2georef -i polygon.geojson -r 4 # vector2georef -r <resolution[0..10]>
> vector2mgrs -i polygon.geojson -r 3 # vector2mgrs -r <resolution[0..5]>
> vector2tilecode -i polygon.geojson -r 18 -p intersect # vector2tilecode -r <resolution[0..29]>
> vector2quadkey -i polygon.geojson -r 18 -p intersect # vector2quadkey -r <resolution[0..29]>
> vector2maidenhead -i polygon.geojson -r 4 # vector2maidenhead -r <resolution[1..4]>
> vector2gars -i polygon.geojson -r 1 # vector2gars -r <resolution[1..4]>
> vector2digipin -i polygon.geojson -r 10 # vector2digipin -r <resolution>
```
- **Uncompact:**
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/vertor2dggs_uncompact.png">
</div>

- **Compact:**
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/vertor2dggs_compact.png">
</div>

### DGGS Compact

``` bash
> h3compact -i h3.geojson -cellid h3 -fix split # h3compact -i <input> -cellid [optional] -f [output_format] -fix [antimeridian method]
> s2compact -i s2.geojson -cellid s2 -fix split # s2compact -i <input> -cellid [optional] -fix [antimeridian method]
> a5compact -i a5.geojson -cellid a5 -split # a5compact -i <input> -cellid [optional] -split
> rhealpixcompact -i rhealpix.geojson -cellid rhealpix -fix split # rhealpixcompact -i <input> -cellid [optional] -fix [antimeridian method]
> dggalcompact -i dggal_cells.geojson -dggs gnosis # dggalcompact -i <input> -dggs <dggal_type> [-zoneid field] [-f output_format]
> dggridcompact -i dggrid_cells.geojson -dggs ISEA4T -r 10 -split -aggregate # Linux API: dggridcompact -i <input> -dggs <DGGS Type> -r <resolution> [-cellid] [-split] [-aggregate]
> isea4tcompact -i isea4t.geojson -cellid isea4t -fix split # isea4tcompact -i <input> -cellid [optional] -fix [antimeridian method]
> isea3hcompact -i isea3h.geojson -cellid isea3h -fix split # isea3hcompact -i <input> -cellid [optional] -fix [antimeridian method]
> easecompact -i ease.geojson -cellid ease # easecompact -i <input> -cellid [optional]
> qtmcompact -i qtm.geojson -cellid qtm # qtmcompact -i <input> -cellid [optional]
> olccompact -i olc.geojson -cellid olc # olccompact -i <input> -cellid [optional]
> geohashcompact -i geohash.geojson -cellid geohash # geohashcompact -i <input> -cellid [optional]
> tilecodecompact -i tilecode.geojson -cellid tilecode # tilecodecompact -i <input> -cellid [optional]
> quadkeycompact -i quadkey.geojson -cellid quadkey # quadkeycompact -i <input> -cellid [optional]
> digipincompact -i digipin.geojson -cellid digipin # digipincompact -i <input> -cellid [optional]
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggscompact_isea4t.png">
</div>


### DGGS Expand

``` bash
> h3expand -i h3_11.geojson -r 12 -cellid h3 -fix split # h3expand -i <input> -r <higher resolution> -cellid [optional] -fix [antimeridian method]
> s2expand -i s2_18.geojson -r 19 -cellid s2 -fix split # s2expand -i <input> -r <higher resolution> -cellid [optional] -fix [antimeridian method]
> a5expand -i a5_18.geojson -r 19 -cellid a5 -split # a5expand -i <input> -r <higher resolution> -cellid [optional] -split
> rhealpixexpand -i rhealpix_10.geojson -r 11 -cellid rhealpix -fix split # rhealpixexpand -i <input> -r <higher resolution> -cellid [optional] -fix [antimeridian method]
> dggalexpand -i dggal_cells.geojson -dggs gnosis -r 6 # dggalexpand -i <input> -dggs <dggal_type> -r <higher resolution>
> dggridexpand -i dggrid_cells.geojson -dggs ISEA4T -r 11 -split -aggregate # Linux API: dggridexpand -i <input> -dggs <DGGS Type> -r <higher resolution> [-cellid] [-split] [-aggregate]
> isea4texpand -i isea4t_18.geojson -r 19 -cellid isea4t -fix split # isea4texpand -i <input> -r <higher resolution> -cellid [optional] -fix [antimeridian method]
> isea3hexpand -i isea3h_18.geojson -r 19 -cellid isea3h -fix split # isea3hexpand -i <input> -r <higher resolution> -cellid [optional] -fix [antimeridian method]
> easeexpand -i ease_4.geojson -r 5 -cellid ease # easeexpand -i <input> -r <higher resolution> -cellid [optional]
> qtmexpand -i qtm_18.geojson -r 19 -cellid qtm # qtmexpand -i <input> -r <higher resolution> -cellid [optional]
> olcexpand -i olc_8.geojson -r 10 -cellid olc # olcexpand -i <input> -r <higher resolution> -cellid [optional]
> geohashexpand -i geohash_7.geojson -r 8 -cellid geohash # geohashexpand -i <input> -r <higher resolution> -cellid [optional]
> tilecodeexpand -i tilecode_18.geojson -r 19 -cellid tilecode # tilecodeexpand -i <input> -r <higher resolution> -cellid [optional]
> quadkeyexpand -i quadkey_18.geojson -r 19 -cellid quadkey # quadkeyexpand -i <input> -r <higher resolution> -cellid [optional]
> digipinexpand -i digipin.geojson -r 11 -cellid digipin # digipinexpand -i <input> -r <higher resolution> -cellid [optional]
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsexpand_isea4t.png">
</div>

### DGGS Resample

Resample a source DGGS layer to another DGGS type (or resolution): build a target grid over the source footprint, then optionally transfer a numeric attribute by **area-weighted overlap** (default) or **nearest-neighbour** assignment from source cells. Omit `-r` or pass `-1` to pick the target resolution that best matches mean source cell area.

``` bash
> dggsresample -i h3_cells.geojson --from_dggs h3 --to_dggs s2 -r 15 -resample_col elevation -m area_weighted -f geojson # dggsresample -i <source layer> --from_dggs <source type> --to_dggs <target type> [-r <resolution> | -1] [-dggs_col <id column>] [-resample_col <numeric field>] [-m area_weighted|nearest] [-f output_format] [-o output_name] [-fix antimeridian method] [-split] [-aggregate] [--dggrid_options JSON] [--a5_options JSON]
```

<div align="center">
  <img src="https://raw.githubusercontent.com/opengeoshub/vgridtools/main/images/readme/dggsresampling_h32s2.png">
</div>

### DGGS Binning
Binning point layers to DGGS cells with spatial joins and aggregation.

Common options for all `*bin` CLIs: `-i` / `--input` (vector file or URL), `-r` / `--resolution`, `-stats` / `--statistics` (`count`, `min`, `max`, `sum`, `mean`, `median`, `std`, `var`, `range`, `minority`, `majority`, `variety`), `-numeric_col` / `--numeric_col` (required when `stats` ≠ `count`), `-category` / `--category` (optional grouping field), `-f` / `--output_format` (default `gpd`).

``` bash
> h3bin -i point.geojson -r 8 -stats count -numeric_col value -category group -fix split # h3bin -i <input> -r <resolution[0..15]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format] -fix [antimeridian method]
> s2bin -i point.geojson -r 13 -stats count -numeric_col value -category group -fix split # s2bin -i <input> -r <resolution[0..30]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format] -fix [antimeridian method]
> a5bin -i point.geojson -r 18 -stats count -numeric_col value -category group -split # a5bin -i <input> -r <resolution[0..29]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format] -split [-options JSON]
> rhealpixbin -i point.geojson -r 8 -stats count -numeric_col value -category group -fix split # rhealpixbin -i <input> -r <resolution[0..15]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format] -fix [antimeridian method]
> dggalbin -i point.geojson -t gnosis -r 5 -stats count -numeric_col value -category group -split # dggalbin -i <input> -t <dggal_type> -r <resolution> -stats [...] -numeric_col [optional] -category [optional] -f [output_format] -split
# dggridbin: Linux API — same options as dggalbin plus -aggregate
> isea4tbin -i point.geojson -r 13 -stats count -numeric_col value -category group -fix split # isea4tbin -i <input> -r <resolution[0..25]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format] -fix [antimeridian method]
> easebin -i point.geojson -r 4 -stats count -numeric_col value -category group # easebin -i <input> -r <resolution[0..6]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format]
> qtmbin -i point.geojson -r 13 -stats count -numeric_col value -category group # qtmbin -i <input> -r <resolution[1..24]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format]
> olcbin -i point.geojson -r 9 -stats count -numeric_col value -category group # olcbin -i <input> -r <resolution[2,4,6,8,10..15]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format]
> geohashbin -i point.geojson -r 6 -stats count -numeric_col value -category group # geohashbin -i <input> -r <resolution[1..10]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format]
> georefbin -i point.geojson -r 4 -stats count # georefbin -i <input> -r <resolution[0..10]>
> tilecodebin -i point.geojson -r 15 -stats count -numeric_col value -category group # tilecodebin -i <input> -r <resolution[0..29]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format]
> quadkeybin -i point.geojson -r 13 -stats count -numeric_col value -category group # quadkeybin -i <input> -r <resolution[0..29]> -stats [...] -numeric_col [optional] -category [optional] -f [output_format]
> maidenheadbin -i point.geojson -r 4 -stats count # maidenheadbin -i <input> -r <resolution[1..4]>
> garsbin -i point.geojson -r 1 -stats count # garsbin -i <input> -r <resolution[1..4]>
> digipinbin -i point.geojson -r 10 -stats count # digipinbin -i <input> -r <resolution> (Python API; no console script yet)
> polygonbin -i points.geojson -p polygons.geojson -stats count # polygonbin -i <points> -p <polygons> [-stats] [-numeric_col] [-category]
```

DGGRID binning (`dggrid_bin` / `dggridbin` Python API, Linux): `-split`, `-aggregate`, same as `vector2dggrid` and `dggridgen`.

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsbinning_h3.png">
</div>

### Raster to DGGS
Convert raster layers in geographic CRS to DGGS. Omit `-r` to use the resolution nearest to the raster cell size. DGGRID/DGGAL support `-aggregate` for binning mode.

``` bash
> raster2h3 -raster raster.tif # raster2h3 -raster <raster> [-r resolution[0..15]] [-stats ...] [-fix antimeridian method]
> raster2s2 -raster raster.tif # raster2s2 -raster <raster> [-r resolution[0..30]] [-fix antimeridian method]
> raster2a5 -raster raster.tif # raster2a5 -raster <raster> [-r resolution[0..29]] [-split]
> raster2rhealpix -raster raster.tif # raster2rhealpix -raster <raster> [-r resolution[0..15]] [-fix antimeridian method]
> raster2dggal -raster raster.tif -dggs gnosis # raster2dggal -raster <raster> -dggs <dggal_type> [-r resolution]
> raster2dggrid -raster raster.tif -t ISEA4T # Linux: raster2dggrid -raster <raster> -t <DGGS Type> [-r resolution] [-split] [-aggregate]
> raster2isea4t -raster raster.tif # raster2isea4t -raster <raster> [-r resolution] (Windows)
> raster2qtm -raster raster.tif # raster2qtm -raster <raster> [-r resolution[1..24]]
> raster2olc -raster raster.tif # raster2olc -raster <raster> [-r resolution]
> raster2geohash -raster raster.tif # raster2geohash -raster <raster> [-r resolution[1..10]]
> raster2georef -raster raster.tif # raster2georef -raster <raster> [-r resolution[0..10]]
> raster2tilecode -raster raster.tif # raster2tilecode -raster <raster> [-r resolution[0..29]]
> raster2quadkey -raster raster.tif # raster2quadkey -raster <raster> [-r resolution[0..29]]
> raster2maidenhead -raster raster.tif # raster2maidenhead -raster <raster> [-r resolution[1..4]]
> raster2gars -raster raster.tif # raster2gars -raster <raster> [-r resolution[1..4]]
> raster2digipin -raster raster.tif # raster2digipin -raster <raster> [-r resolution]
```

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/raster2dggs_h3.png">
</div>


## DGGS Generator
Generate DGGS at a specific resolution within a bounding box (`-b min_lon min_lat max_lon max_lat`). Antimeridian: `-fix` (H3, S2, rHEALPix, ISEA4T/ISEA3H), `-split` (A5, DGGAL, DGGRID), `-aggregate` (DGGRID).

``` bash
> h3grid -r 11 -b 106.699007 10.762811 106.717674 10.778649 -fix split # h3grid -r <resolution[0..15]> -b <bbox> [-fix method] [-f output_format]
> s2grid -r 18 -b 106.699007 10.762811 106.717674 10.778649 -fix split # s2grid -r <resolution[0..30]> -b <bbox> [-fix method]
> a5grid -r 18 -b 106.699007 10.762811 106.717674 10.778649 -split # a5grid -r <resolution[0..29]> -b <bbox> [-split] [-options JSON]
> rhealpixgrid -r 11 -b 106.699007 10.762811 106.717674 10.778649 -fix split # rhealpixgrid -r <resolution[0..15]> -b <bbox> [-fix method]
> dggalgen -dggs gnosis -r 5 -bbox 106.699007,10.762811,106.717674,10.778649 # dggalgen -dggs <dggal_type> -r <resolution> [-bbox min_lon,min_lat,max_lon,max_lat] [-c compact]
> dggridgen -t ISEA3H -r 2 -b 106.699007 10.762811 106.717674 10.778649 -split -aggregate # Linux: dggridgen -t <DGGS Type> -r <resolution> [-b bbox] [-a address_type] [-split] [-aggregate]
> isea4tgrid -r 17 -b 106.699007 10.762811 106.717674 10.778649 -fix split # isea4tgrid -r <resolution[0..25]> -b <bbox> [-fix method] (Windows)
> isea3hgrid -r 20 -b 106.699007 10.762811 106.717674 10.778649 -fix split # isea3hgrid -r <resolution[0..40]> -b <bbox> [-fix method] (Windows)
> easegrid -r 4 -b 106.699007 10.762811 106.717674 10.778649 # easegrid -r <resolution[0..6]> -b <bbox>
> qtmgrid -r 8 -b 106.699007 10.762811 106.717674 10.778649 # qtmgrid -r <resolution[1..24]> -b <bbox>
> olcgrid -r 8 -b 106.699007 10.762811 106.717674 10.778649 # olcgrid -r <resolution[2,4,6,8,10..15]> -b <bbox>
> geohashgrid -r 6 -b 106.699007 10.762811 106.717674 10.778649 # geohashgrid -r <resolution[1..10]> -b <bbox>
> georefgrid -r 4 -b 106.699007 10.762811 106.717674 10.778649 # georefgrid -r <resolution[0..10]> -b <bbox>
> mgrsgrid -r 1 -gzd 48P # mgrsgrid -r <resolution[0..5]> -gzd <Grid Zone Designator, e.g. 48P>
> tilecodegrid -r 20 -b 106.699007 10.762811 106.717674 10.778649 # tilecodegrid -r <resolution[0..29]> -b <bbox>
> quadkeygrid -r 20 -b 106.699007 10.762811 106.717674 10.778649 # quadkeygrid -r <resolution[0..29]> -b <bbox>
> maidenheadgrid -r 4 -b 106.699007 10.762811 106.717674 10.778649 # maidenheadgrid -r <resolution[1..4]> -b <bbox>
> garsgrid -r 1 -b 106.699007 10.762811 106.717674 10.778649 # garsgrid -r <resolution[1..4]> -b <bbox>
> digipingrid -r 10 -b 106.699007 10.762811 106.717674 10.778649 # digipingrid -r <resolution> -b <bbox>
```

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsgenerator_h3.png">
</div>

## DGGS Inspect
``` bash
> h3inspect -r 3 # h3inspect -r <resolution>
> s2inspect -r 6
> a5inspect -r 5
> rhealpixinspect -r 4
> dggalinspect -t gnosis -r 5 # dggalinspect -t <dggal_type> -r <resolution>
> dggridinspect -t ISEA4T -r 3 # Linux: dggridinspect -t <DGGS Type> -r <resolution>
> isea4tinspect -r 5
> isea3hinspect -r 3
> easeinspect -r 0
> qtminspect -r 6
> olcinspect -r 4
> geohashinspect -r 3
> georefinspect -r 1
> tilecodeinspect -r 7
> quadkeyinspect -r 7
> maidenheadinspect -r 2
> garsinspect -r 1
> digipininspect -r 10
```

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/a5_inspect.png">
</div>

### Distribution of DGGS Area Distortions visualized from DGGS Inspect
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/a5_norm_area.png">
</div>


### Visualization of DGGS IPQ Compactness visualized from DGGS Inspect
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/a5_compactness.png">
</div>


## DGGS Stats

``` bash
> h3stats # Number of cells, avg edge length, avg cell area at each resolution
> s2stats
> a5stats
> rhealpixstats
> dggalstats -t gnosis # dggalstats -t <dggal_type>
> dggridstats -t FULLER3H -r 8 # Linux: dggridstats -t <DGGS Type> [-r resolution] [-aggregate]
> isea4tstats
> isea3hstats
> easestats
> qtmstats
> olcstats
> geohashstats
> georefstats
> mgrsstats
> tilecodestats
> quadkeystats
> maidenheadstats
> garsstats
> digipinstats
```

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsstats.png">
</div>


## Polyhedra Generator
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/polyhedra.png">
</div>

``` bash
> tetrahedron  # Generate global tetrahedron
> cube         # Generate global cube
> octahedron   # Generate global octahedron
> dodecahedron # Generate global dodecahedron
> fuller_icosahedron   # Generate global Fuller icosahedron
> rhombic_icosahedron   # Generate global rhombic icosahedron (DGGAL)
``` 

# Vgrid
**Vgrid -  A unified framework for working with DGGS**

Vgrid supports a wide range of popular DGGSs, including H3, S2, A5, rHEALPix, DGGAL, DGGRID, Open-EAGGR ISEA4T, ISEA3H, EASE-DGGS, QTM, OLC, Geohash, GEOREF, MGRS, TileCode, Quadkey, Maidenhead and GARS.

Vgrid supports popular input and output GIS formats, including CSV, GeoJSON, Pandas/GeoPandas, Shapefile, GeoPackage, and GeoParquet.

[![logo](https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs.png)](https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs.png)

Full Vgrid DGGS documentation is available at [vgrid document](https://vgrid.gishub.vn).

To work with Vgrid DGGS directly in GeoPandas and Pandas, use the [vgridpandas](https://pypi.org/project/vgridpandas/) package. Full Vgridpandas DGGS documentation is available at [vgridpandas document](https://vgridpandas.gishub.vn).

To work with Vgrid DGGS in QGIS, install the [Vgrid Plugin](https://plugins.qgis.org/plugins/vgridtools/).

To visualize DGGS in Maplibre GL JS, try the [vgrid-maplibre](https://www.npmjs.com/package/vgrid-maplibre) library.

For an interactive demo, visit the [Vgrid Homepage](https://vgrid.vn).

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opengeoshub/vgrid/blob/main)
[![image](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/opengeoshub/vgrid/HEAD?filepath=docs/notebooks)
[![image](https://studiolab.sagemaker.aws/studiolab.svg)](https://studiolab.sagemaker.aws/import/github/opengeoshub/vgrid/blob/main/docs/notebooks/00_intro.ipynb)
[![image](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://jupyterlite.gishub.vn/lab/index.html?path=notebooks/vgrid/00_intro.ipynb)
[![PyPI version](https://badge.fury.io/py/vgrid.svg)](https://badge.fury.io/py/vgrid)
[![image](https://static.pepy.tech/badge/vgrid)](https://pepy.tech/project/vgrid)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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

##  Area distortion over normalized areas of popular geodesic DGGS
<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_norm_area_box.png">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_norm_area.png">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_norm_area_distribution.png">
</p>


##  IPQ compactness distribution over popular geodesic DGGS

Isoperimetric Inequality (IPQ) Compactness (suggested by [Osserman, 1978](https://sites.math.washington.edu/~toro/Courses/20-21/MSF/osserman.pdf)):

$$C_{IPQ} = \frac{4 \pi A}{p^2}$$

The range of the IPQ compactness metric is [0,1]. 

A circle represents the maximum compactness with a value of 1. 

As shapes become more irregular or elongated, their compactness decreases toward 0.

<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_ipq_box.png">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_ipq.png">
</p>

##  Convex hull compactness distribution over popular geodesic DGGS

$$C_{CVH} = \frac{A}{A_{CVH}}$$

The range of the convex hull compactness metric is [0,1]. 

As shapes become more concave, their convex hull compactness decreases toward 0.
<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_cvh_box.png">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_cvh.png">
</p>

## Installation
### pip
[![image](https://img.shields.io/pypi/v/vgrid.svg)](https://pypi.python.org/pypi/vgrid)
```bash
pip install vgrid --upgrade
```

## Key Features

- **DGGS Conversion:** Convert Latlon to DGGS, DGGS cells to Shapely Geometry/ GeoJSON, Vector to DGGS, Raster to DGGS.

- **DGGS Compact:** Compact DGGS cells or expand them to a specific resolution.

- **DGGS Resample:** Resample a source DGGS layer to another DGGS type or resolution, with optional area-weighted or nearest-neighbour attribute transfer.

- **DGGS Binning:** Aggregate points into DGGS cells, supporting common statistics (count, min, max, etc.) and category-based groups.

- **DGGS Generator:** Generate DGGS at a specfic bounding box and resolution.

- **DGGS Inspect:** Calculate and visualize DGGS area distortions and IPQ (isoperimetric inequality) compactness at a specific resolution.

- **DGGS Stats:** Show DGGS metrics for each resolution like number of cells, average edge length, average cell area, perimeter.

## Usage examples

### Latlon to DGGS

```python
>>> from vgrid.conversion.latlon2dggs import latlon2h3
>>> lat = 10.775276
>>> lon = 106.706797
>>> res = 10
h3_id = latlon2h3(lat, lon, res)
h3_id

'8a65b56628e7fff'
```

### DGGS to Shapely Polygon

```python
import geopandas as gpd
from vgrid.conversion.dggs2geo.h32geo import h32geo
h3_geo = h32geo(h3_id)
h3_geo

POLYGON ((106.70713615426936 10.774978441229653, 106.70721514572995 10.775713374905791, 106.70661718075799 10.776150587194028, 106.70594022237229 10.77585286364977, 106.70586123315323 10.77511792678204, 106.70645920007833 10.774680716650128, 106.70713615426936 10.774978441229653))
```

### DGGS to GeoJSON

```python
>>> from vgrid.conversion.dggs2geo.h32geo import h32geojson
>>> h3_geojson = h32geojson(h3_id)
>>> h3_geojson

{'type': 'FeatureCollection', 'features': [{'type': 'Feature', 'geometry': {'type': 'Polygon', 'coordinates': (((106.70713615426936, 10.774978441229653), (106.70721514572995, 10.775713374905791), (106.70661718075799, 10.776150587194028), (106.70594022237229, 10.77585286364977), (106.70586123315323, 10.77511792678204), (106.70645920007833, 10.774680716650128), (106.70713615426936, 10.774978441229653)),)}, 'properties': {'h3': '8a65b56628e7fff', 'resolution': 10, 'center_lat': 10.7754157, 'center_lon': 106.7065382, 'avg_edge_len': 81.374, 'cell_area': 17202.984}}]}
```

### Vector to DGGS

```python
from vgrid.conversion.vector2dggs.vector2isea4t import vector2isea4t

file_path = ("https://raw.githubusercontent.com/opengeoshub/vopendata/main/shape/polygon.geojson")
vector_to_isea4t = vector2isea4t(file_path, resolution=16, compact=False, predicate = "centroid_within", output_format="gpd")
```

### DGGS Compact

```python
from vgrid.conversion.dggscompact.isea4tcompact import isea4tcompact

isea4t_compacted = isea4tcompact(vector_to_isea4t,output_format="gpd")
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggscompact_isea4t.png">
</div>


### DGGS Expand

```python
from vgrid.conversion.dggscompact.isea4tcompact import isea4texpand

isea4t_expanded = isea4texpand(isea4t_compacted, resolution=17, output_format="gpd")   
```

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsexpand_isea4t.png">
</div>

### DGGS Resample

Resample a source DGGS layer to another DGGS type (or resolution): build a target grid over the source footprint, then optionally transfer a numeric attribute by **area-weighted overlap** (default) or **nearest-neighbour** assignment from source cells. Omit `resolution` or pass `-1` to pick the target resolution that best matches mean source cell area.

```python
from vgrid.conversion.dggsresample.dggsresample import dggsresample
from vgrid.conversion.vector2dggs.vector2h3 import vector2h3

file_path = "https://raw.githubusercontent.com/opengeoshub/vopendata/main/shape/polygon.geojson"
h3_cells = vector2h3(file_path, resolution=10, output_format="gpd")
s2_resampled = dggsresample(
    h3_cells,
    dggs_from="h3",
    dggs_to="s2",
    resolution=15,
    method="area_weighted",
    output_format="gpd",
)
```

<div align="center">
  <img src="https://raw.githubusercontent.com/opengeoshub/vgridtools/main/images/readme/dggsresampling_h32s2.png">
</div>

### DGGS Binning

```python
from vgrid.binning.h3bin import h3bin
file_path = ("https://raw.githubusercontent.com/opengeoshub/vopendata/main/csv/dist1_pois.csv")
stats="count"
h3_bin = h3bin(file_path, resolution=10, stats=stats, 
                # numeric_col="confidence",
                # category="category",
                output_format="gpd")
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsbinning_h3.png">
</div>


### Raster to DGGS

```python
from vgrid.conversion.raster2dggs.raster2h3 import raster2h3
from vgrid.utils.io import download_file

raster_url = ("https://raw.githubusercontent.com/opengeoshub/vopendata/main/raster/rgb.tif")
raster_file = download_file(raster_url)
raster_to_h3 =  raster2h3(raster_file,output_format="gpd")
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/raster2dggs_h3.png">
</div>


### DGGS Generator

```python
from vgrid.generator.h3grid import h3grid

h3_grid = h3grid(resolution=0, output_format="gpd")
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsgenerator_h3.png">
</div>


## DGGS Inspect
```python
from vgrid.stats.a5stats import a5inspect
a5inspect()
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

### DGGS Stats
```python
from vgrid.stats.h3stats import h3stats
h3stats()
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsstats.png">
</div>

### Further examples
For more examples, see the 
[example notebooks](https://nbviewer.jupyter.org/github/opengeoshub/vgrid/tree/main/docs/notebooks/).
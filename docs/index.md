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
[![image](https://studiolab.sagemaker.aws/studiolab.svg)](https://studiolab.sagemaker.aws/import/github/opengeoshub/vgrid/blob/main/docs/notebooks/01_h3.ipynb)
[![image](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://jupyterlite.gishub.vn/lab/index.html?path=notebooks/vgrid/01_h3.ipynb)
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

The range of the IPQ compactness metric is (0,1]. 

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

The range of the convex hull compactness metric is (0,1]. 

As shapes become more concave, their convex hull compactness decreases toward 0.
<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_cvh_box.png">
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs_cvh.png">
</p>
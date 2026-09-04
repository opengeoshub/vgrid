# Vgrid DGGS key features

- **DGGS Conversion:** Convert Latlon to DGGS, DGGS cells to Shapely Geometry/ GeoJSON, Vector to DGGS, Raster to DGGS.

- **DGGS Compact:** Compact DGGS cells or expand them to a specific resolution.

- **DGGS Resample:** Resample a source DGGS layer to another DGGS type or resolution, with optional area-weighted or nearest-neighbour attribute transfer and a source–target keep filter (`centroid_within` or `intersects`).

- **DGGS Binning:** Aggregate points into DGGS cells, supporting common statistics (count, min, max, etc.) and category-based groups.

- **DGGS Generator:** Generate DGGS at a specfic bounding box and resolution.

- **DGGS Inspect:** Calculate and visualize DGGS area distortions and IPQ (isoperimetric inequality) compactness at a specific resolution.

- **DGGS Stats:** Show DGGS metrics for each resolution like number of cells, average edge length, average cell area, perimeter.

## Usage examples

### Latlon to DGGS

```python
from vgrid.conversion.latlon2dggs import latlon2h3
lat = 10.775276
lon = 106.706797
res = 10
h3_id = latlon2h3(lat, lon, res)
```

### DGGS to Shapely Polygon

```python
import geopandas as gpd
from vgrid.conversion.dggs2geo.h32geo import h32geo
h3_geo = h32geo(h3_id)
```

### DGGS to GeoJSON

```python
from vgrid.conversion.dggs2geo.h32geo import h32geojson
h3_geojson = h32geojson(h3_id)
```

### Vector to DGGS

```python
from vgrid.conversion.vector2dggs.vector2isea4t import vector2isea4t

file_path = ("https://raw.githubusercontent.com/opengeoshub/vopendata/main/shape/polygon.geojson")

vector_to_isea4t = vector2isea4t(
    file_path,
    resolution=16,
    compact=False,
    depth=-1,  # used when compact=True; -1: full compact
    predicate="centroid_within",
    output_format="gpd",
    verbose=True,
)
```

### DGGS Compact

```python
from vgrid.conversion.dggscompact.isea4tcompact import isea4tcompact

isea4t_compacted = isea4tcompact(
    vector_to_isea4t,
    depth=-1,  # -1: full compact
    output_format="gpd",
    verbose=True,
)
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggscompact_isea4t.png">
</div>


### DGGS Expand

```python
from vgrid.conversion.dggscompact.isea4tcompact import isea4texpand

isea4t_expanded = isea4texpand(
    isea4t_compacted,
    resolution=17,  # if set, depth is ignored
    # depth=1,  # 1: children
    output_format="gpd",
    verbose=True,
)
```

<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsexpand_isea4t.png">
</div>

### DGGS Resample

Resample a source DGGS layer to another DGGS type (or resolution): build a target grid over the source footprint, then optionally transfer a numeric attribute by **area-weighted overlap** (default) or **nearest-neighbour** assignment from source cells. Keep target cells with **`centroid_within`** (default: the target cell contains a source centroid) or **`intersects`**. Omit `resolution` or pass `-1` to pick the target resolution that best matches mean source cell area.

```python
from vgrid.conversion.dggsresample.dggsresample import dggsresample
from vgrid.conversion.vector2dggs.vector2h3 import vector2h3

file_path = "https://raw.githubusercontent.com/opengeoshub/vopendata/main/shape/polygon.geojson"
h3_cells = vector2h3(file_path, resolution=10, output_format="gpd", verbose=True)
s2_resampled = dggsresample(
    h3_cells,
    dggs_from="h3",
    dggs_to="s2",
    resolution=15,
    method="area_weighted",
    predicate="centroid_within",
    output_format="gpd",
    verbose=True,
)
```

<div align="center">
  <img src="https://raw.githubusercontent.com/opengeoshub/vgridtools/main/images/readme/dggsresampling_h32s2.png">
</div>

### DGGS Binning

```python
from vgrid.binning.h3bin import h3bin
file_path = ("https://raw.githubusercontent.com/opengeoshub/vopendata/main/csv/dist1_pois.csv")
agg="count"
h3_bin = h3bin(file_path, resolution=10, agg=agg, 
                # numeric_col="confidence",
                # category="category",
                output_format="gpd", verbose=True)
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
raster_to_h3 = raster2h3(raster_file, output_format="gpd", verbose=True)
```
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/raster2dggs_h3.png">
</div>


### DGGS Generator

```python
from vgrid.generator.h3grid import h3grid

h3_grid = h3grid(resolution=0, output_format="gpd", verbose=True)
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


### Distribution of DGGS IPQ Compactness visualized from DGGS Inspect
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

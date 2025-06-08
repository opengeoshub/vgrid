# Usage

You can try out vgrid by using Goolge Colab ([![image](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opengeoshub/vgrid/blob/master)) or Binder ([![image](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/opengeoshub/vgrid/HEAD)) without having to install anything on your computer.

## Launch Jupyter notebook

```bash
conda activate geo
jupyter notebook
```

## Use ipyleaflet plotting backend

```python
import vgrid
```

## Use folium plotting backend

```python
import vgrid.foliumap as vgrid
```

## Create an interactive map

```python
m = vgrid.Map(center=(40, -100), zoom=4)
m
```

## Customize map height

```python
m = vgrid.Map(height="450px")
m
```

## Set control visibility

```python
m = vgrid.Map(draw_control=False, measure_control=False, fullscreen_control=False, attribution_control=True)
m
```

## Change basemaps

```python
m = vgrid.Map(google_map="TERRAIN")
m.add_basemap("HYBRID")
m
```

## Add XYZ tile layer

```python
m = vgrid.Map()
m.add_tile_layer(url="https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}", name="Google Satellite", attribution="Google")
m
```

## Add WMS tile layer

```python
m = vgrid.Map()
naip_url = 'https://services.nationalmap.gov/arcgis/services/USGSNAIPImagery/ImageServer/WMSServer?'
m.add_wms_layer(url=naip_url, layers='0', name='NAIP Imagery', format='image/png', shown=True)
m
```

## Use HERE Map Widget for Jupyter plotting backend

### Prerequisites

-   A HERE developer account, free and available under [HERE Developer Portal](https://developer.here.com)
-   An [API key](https://developer.here.com/documentation/identity-access-management/dev_guide/topics/dev-apikey.html) from the [HERE Developer Portal](https://developer.here.com)
-   Export API key into environment variable `HEREMAPS_API_KEY`

```bash
export HEREMAPS_API_KEY=YOUR-ACTUAL-API-KEY
```

```python
import vgrid.heremap as vgrid
```

### Create an interactive map

```python
import os
api_key = os.environ.get("HEREMAPS_API_KEY") # read api_key from environment variable.
m = vgrid.Map(api_key=api_key, center=(40, -100), zoom=4)
m
```

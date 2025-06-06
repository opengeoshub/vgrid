# Vgrid Tutorial

This tutorial demonstrates the main features and functionality of the `vgrid` package, a comprehensive toolkit for working with Discrete Global Grid Systems (DGGS) and cell-based geocoding.

## Installation

First, let's install the package:

```bash
pip install vgrid --upgrade
```

---

## 1. DGGS Conversion

### 1.1 Converting Latitude/Longitude to DGGS

Let's start with converting coordinates to different DGGS formats:

```python
from vgrid.conversion import latlon2h3, latlon2s2, latlon2rhealpix

# Example coordinates (Ho Chi Minh City)
lat, lon = 10.775276, 106.706797

# Convert to H3
h3_cell = latlon2h3(lat, lon, 13)
print(f"H3 cell: {h3_cell}")

# Convert to S2
s2_cell = latlon2s2(lat, lon, 21)
print(f"S2 cell: {s2_cell}")

# Convert to rHEALPix
rhealpix_cell = latlon2rhealpix(lat, lon, 14)
print(f"rHEALPix cell: {rhealpix_cell}")
```

---

### 1.2 Converting DGGS to GeoJSON

Let's convert DGGS cells to GeoJSON format for visualization:

```python
from vgrid.conversion import h32geojson, s22geojson
import json

# Convert H3 cell to GeoJSON
h3_geojson = h32geojson(h3_cell)
print("H3 GeoJSON:")
print(json.dumps(h3_geojson, indent=2))

# Convert S2 cell to GeoJSON
s2_geojson = s22geojson(s2_cell)
print("\nS2 GeoJSON:")
print(json.dumps(s2_geojson, indent=2))
```

---

### 1.3 Vector to DGGS Conversion

Let's convert a polygon to DGGS cells:

```python
from vgrid.conversion import geojson2h3

# Create a sample polygon
polygon = {
    "type": "Feature",
    "geometry": {
        "type": "Polygon",
        "coordinates": [[
            [106.7, 10.7],
            [106.8, 10.7],
            [106.8, 10.8],
            [106.7, 10.8],
            [106.7, 10.7]
        ]]
    }
}

# Convert to H3 cells
h3_cells = geojson2h3(polygon, resolution=11, compact=True)
print(f"Number of H3 cells: {len(h3_cells)}")
```

---

## 2. DGGS Operations

### 2.1 DGGS Compaction

Let's demonstrate how to compact DGGS cells:

```python
from vgrid.conversion import h3compact

# Compact H3 cells
compacted_cells = h3compact(h3_cells)
print(f"Number of cells before compaction: {len(h3_cells)}")
print(f"Number of cells after compaction: {len(compacted_cells)}")
```

---

### 2.2 DGGS Expansion

Let's demonstrate how to expand DGGS cells to a higher resolution:

```python
from vgrid.conversion import h3expand

# Expand H3 cells to a higher resolution
expanded_cells = h3expand(h3_cells, resolution=12)
print(f"Number of cells before expansion: {len(h3_cells)}")
print(f"Number of cells after expansion: {len(expanded_cells)}")
```

---

## 3. DGGS Binning

Let's demonstrate how to bin points into DGGS cells:

```python
from vgrid.binning import h3bin
import pandas as pd
import json

# Create sample points
points = pd.DataFrame({
    'lat': [10.7, 10.75, 10.8],
    'lon': [106.7, 106.75, 106.8],
    'value': [1, 2, 3]
})

# Convert to GeoJSON
points_geojson = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [row['lon'], row['lat']],
            },
            "properties": {"value": row['value']},
        }
        for _, row in points.iterrows()
    ]
}

# Bin points into H3 cells
binned_cells = h3bin(points_geojson, resolution=8, stats='count', field='value')
print("Binned cells:")
print(json.dumps(binned_cells, indent=2))
```

---

## 4. DGGS Resampling

Let's demonstrate how to resample between different DGGS types:

```python
from vgrid.resampling import resample

# Create sample H3 data
h3_data = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "geometry": h3_geojson['geometry'],
            "properties": {"population": 1000}
        }
    ]
}

# Resample from H3 to S2
resampled_data = resample(h3_data, from_dggs='h3', to_dggs='s2', resample_field='population')
print("Resampled data:")
print(json.dumps(resampled_data, indent=2))
```

---

## 5. Visualization

Let's visualize some of our DGGS cells using folium:

```python
import folium
from vgrid.conversion import h32geojson

# Create a map centered on Ho Chi Minh City
m = folium.Map(location=[10.775276, 106.706797], zoom_start=12)

# Add H3 cells
for cell in h3_cells:
    cell_geojson = h32geojson(cell)
    folium.GeoJson(
        cell_geojson,
        style_function=lambda x: {'fillColor': 'blue', 'color': 'black', 'weight': 1, 'fillOpacity': 0.3}
    ).add_to(m)

# Display the map
m
```

---

## Tips
- You may need to install extra dependencies:
  ```bash
  pip install folium geopandas pandas
  ```
- If you want to add more examples (e.g., raster conversion, other DGGS types), just add more code blocks following the same pattern. 
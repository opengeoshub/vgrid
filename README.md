# Vgrid - Global Geocoding Systems

## Installation: 
- Using pip:   
    ``` bash 
    pip install vgrid --upgrade
    ```
    
- Visit vgrid on [PyPI](https://pypi.org/project/vgrid/)

## Usage - Python:
### Import vgrid, initialize Latitude and Longitude for testing:
``` python
from vgrid.geocode import olc, mgrs, geohash, georef, s2, h3, maidenhead, gars 
import h3, json
from vgrid.geocode.gars import GARSGrid
from vgrid.geocode.s2 import LatLng, CellId
from vgrid.geocode.geocode2geojson import *

latitude, longitude = 10.775275567242561, 106.70679737574993
print(f'Latitude, Longitude: ({latitude}, {longitude})')
```

### OLC
``` python
print('\OLC:')
olc_precision = 11 #[10-->15]
olc_code = olc.encode(latitude, longitude, olc_precision)
olc_decode = olc.decode(olc_code)
print(f'OLC at precision = {olc_precision}: {olc_code}')
print(f'Decode {olc_code} to center and cell in WGS84 = {olc_decode}')

data = olc2geojson(olc_code)
output_file = f'olc{olc_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### MGRS
``` python
print('\nMGRS:')
mgrs_precision = 4 # [0 -->5]
mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_precision)
mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
print(f'MGRS Code at precision = {mgrs_precision}: {mgrs_code}')
print(f'Convert {mgrs_code} to WGS84 = {mgrs_code_to_wgs}')

data = mgrs2geojson(mgrs_code)
output_file = f'mgrs{mgrs_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### Geohash
``` python
print('\Geohsash:')
geohash_precision = 9 # [1-->30]
geohash_code = geohash.encode(latitude, longitude, geohash_precision)
geohash_decode = geohash.decode(geohash_code, True)
print(f'Geohash Code at precision = {geohash_precision}: {geohash_code}')
print(f'Decode {geohash_code} to WGS84 = {geohash_decode}')

data = geohash2geojson(geohash_code)
output_file = f'geohash{geohash_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  
print(f'GeoJSON written to {output_file}')
```

### GEOREF
``` python
print('\GEOREF:')
georef_precision = 4 # [0 -->10]
georef_code = georef.encode(latitude, longitude, georef_precision)
georef_decode = georef.decode(georef_code, True)
print(f'latitude, longitude = {latitude},{longitude}')

print(f'GEOREF Code at precision = {georef_precision}: {georef_code}')
print(f'Decode {georef_code} to WGS84 = {georef_decode}')

data = georef2geojson(georef_code)
output_file = f'georef{georef_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### S2
``` python
print('\S2:')
s2_precision = 21 #[0 -->30]
lat_lng = LatLng.from_degrees(latitude, longitude)
cell_id = CellId.from_lat_lng(lat_lng)
cell_id = cell_id.parent(s2_precision)
cell_id_token= CellId.to_token(cell_id)
print(cell_id_token)

data = s22geojson(cell_id_token)
output_file = f's2_{s2_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### H3
``` python
print('\H3:')
h3_precision = 13 #[0-->15]
h3_code = h3.latlng_to_cell(latitude, longitude, h3_precision)
h3_decode = h3.cell_to_latlng(h3_code)

print(f'latitude, longitude = {latitude},{longitude}')
print(f'H3 code at precision = {h3_precision}: {h3_code}')
print(f'Decode {h3_code} to WGS84 = {h3_decode}')

data = h32geojson(h3_code)
output_file = f'h3_{h3_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### Maidenhead
``` python
print('\nMaidenhead:')
maidenhead_precision = 4 #[1-->4]
maidenhead_code = maidenhead.toMaiden(latitude, longitude, maidenhead_precision)
maidenGrid = maidenhead.maidenGrid(maidenhead_code)
print(f'Maidenhead Code at precision = {maidenhead_precision}: {maidenhead_code}')
print(f'Convert {maidenhead_code} to center and cell in WGS84 = {maidenGrid}')

data = maidenhead2geojson(maidenhead_code)
output_file = f'maidenhead_{maidenhead_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### GARS
``` python
print('\nGARS:')
gars_precision = 1 # 1, 5, 15, 30 minutes
gars_grid = GARSGrid.from_latlon(latitude, longitude, gars_precision)
gars_code = gars_grid.gars_id
print(gars_code)

data = gars2geojson(gars_code)
output_file = f'gars_{gars_precision}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  
print(f'GeoJSON written to {output_file}')
```
## Usage - CLI:
### OLC
``` bash
> olc2geojson 7P28QPG4+4P7
```

### MGRS
``` bash
> mgrs2geojson 34TGK56063228
```

### Geohash
``` bash
> geohash2geojson w3gvk1td8
```

### GEOREF
``` bash
> georef2geojson VGBL42404651
```

### S2
``` bash
> s22geojson 31752f45cc94 
```

### H3
``` bash
> h32geojson 8d65b56628e46bf 
```

### Maidenhead
``` bash
> maidenhead2geojson OK30is46 
```

### GARS
``` bash
> gars2geojson 574JK1918
```

### Command line for creating geocoding grids in shapefile format
``` bash
> geohashgrid -p 1 -o geohash_1.shp (p = [1..12])
> maidenheadgrid -p 1 -o maidenhead_1.shp (p = [1, 2, 3, 4])
> gzd -o gzd.shp (Create Grid Zone Designators - used by MGRS)
> mgrsgrid -o mgrs_32648.shp -cellsize 100000 -epsg 32648 (Create MGRS Grid with cell size 100km x 100km at UTM zone 48N)  
```
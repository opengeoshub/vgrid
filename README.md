# Vgrid - DGGS and Cell-based Geocoding Utilities

## Installation: 
- Using pip:   
    ``` bash 
    pip install vgrid --upgrade
    ```
    
- Visit vgrid on [PyPI](https://pypi.org/project/vgrid/)

## Demo Page:  [Vgrid Home](https://vgrid.vn)

## Usage - Python:
### Import vgrid, initialize latitude and longitude for testing:
``` python
from vgrid.geocode import h3,s2, olc, geohash, georef, mgrs, tilecode, maidenhead, gars 
import h3, json
from vgrid.utils.gars.garsgrid import GARSGrid
from vgrid.geocode.s2 import LatLng, CellId
from vgrid.geocode.geocode2geojson import *
from vgrid.geocode.latlon2geocode import *


latitude, longitude = 10.775275567242561, 106.70679737574993
print(f'Latitude, Longitude: ({latitude}, {longitude})')
```

### H3
``` python
print('\H3:')
h3_resolution = 13 #[0..15]
h3_code = h3.latlng_to_cell(latitude, longitude, h3_resolution)
h3_decode = h3.cell_to_latlng(h3_code)

print(f'latitude, longitude = {latitude},{longitude}')
print(f'H3 code at resolution = {h3_resolution}: {h3_code}')
print(f'Decode {h3_code} to WGS84 = {h3_decode}')

data = h32geojson(h3_code)
output_file = f'h3_{h3_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### S2
``` python
print('\S2:')
s2_resolution = 21 #[0..30]
lat_lng = LatLng.from_degrees(latitude, longitude)
cell_id = CellId.from_lat_lng(lat_lng)
cell_id = cell_id.parent(s2_resolution)
cell_id_token= CellId.to_token(cell_id)
print(cell_id_token)

data = s22geojson(cell_id_token)
output_file = f's2_{s2_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```


### OLC
``` python
print('\OLC:')
olc_resolution = 11 #[10..15]
olc_code = olc.encode(latitude, longitude, olc_resolution)
olc_decode = olc.decode(olc_code)
print(f'OLC at resolution = {olc_resolution}: {olc_code}')
print(f'Decode {olc_code} to center and cell in WGS84 = {olc_decode}')

data = olc2geojson(olc_code)
output_file = f'olc{olc_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### Geohash
``` python
print('\Geohsash:')
geohash_resolution = 9 # [1..30]
geohash_code = geohash.encode(latitude, longitude, geohash_resolution)
geohash_decode = geohash.decode(geohash_code, True)
print(f'Geohash Code at resolution = {geohash_resolution}: {geohash_code}')
print(f'Decode {geohash_code} to WGS84 = {geohash_decode}')

data = geohash2geojson(geohash_code)
output_file = f'geohash{geohash_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  
print(f'GeoJSON written to {output_file}')
```

### GEOREF
``` python
print('\GEOREF:')
georef_resolution = 4 # [0..10]
georef_code = georef.encode(latitude, longitude, georef_resolution)
georef_decode = georef.decode(georef_code, True)
print(f'latitude, longitude = {latitude},{longitude}')

print(f'GEOREF Code at resolution = {georef_resolution}: {georef_code}')
print(f'Decode {georef_code} to WGS84 = {georef_decode}')

data = georef2geojson(georef_code)
output_file = f'georef{georef_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### MGRS
``` python
print('\nMGRS:')
mgrs_resolution = 4 # [0..5]
mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_resolution)
mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
print(f'MGRS Code at resolution = {mgrs_resolution}: {mgrs_code}')
print(f'Convert {mgrs_code} to WGS84 = {mgrs_code_to_wgs}')

data = mgrs2geojson(mgrs_code)
output_file = f'mgrs{mgrs_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### Tilecode
``` python
print('\Tilecode:')
tilecode_esolution = 23  # [0..26]
tile_code = tilecode.latlon2tilecode(latitude, longitude, tilecode_esolution)
tile_encode = tilecode.tilecode2latlon(tile_code)
print(f'Tilecode at zoom level = {resolution}: {tile_code}')
print(f'Convert {tile_code} to WGS84 = {tile_encode}')

data = tilecode2geojson(tile_code)
print(data)

output_file = f'tilecode_{resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```


### Maidenhead
``` python
print('\nMaidenhead:')
maidenhead_resolution = 4 #[1..4]
maidenhead_code = maidenhead.toMaiden(latitude, longitude, maidenhead_resolution)
maidenGrid = maidenhead.maidenGrid(maidenhead_code)
print(f'Maidenhead Code at resolution = {maidenhead_resolution}: {maidenhead_code}')
print(f'Convert {maidenhead_code} to center and cell in WGS84 = {maidenGrid}')

data = maidenhead2geojson(maidenhead_code)
output_file = f'maidenhead_{maidenhead_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  # 'indent' makes the JSON output more readable
print(f'GeoJSON written to {output_file}')
```

### GARS
``` python
print('\nGARS:')
gars_resolution = 1 # [1, 5, 15, 30 minutes]
gars_grid = GARSGrid.from_latlon(latitude, longitude, gars_resolution)
gars_code = gars_grid.gars_id
print(gars_code)

data = gars2geojson(gars_code)
output_file = f'gars_{gars_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  
print(f'GeoJSON written to {output_file}')
```
## Usage - CLI:
### H3
``` bash
> h32geojson 8d65b56628e46bf 
> h3stats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### S2
``` bash
> s22geojson 31752f45cc94 
> s2stats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### OLC
``` bash
> olc2geojson 7P28QPG4+4P7
> olcstats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### Geohash
``` bash
> geohash2geojson w3gvk1td8
> geohashstats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### GEOREF
``` bash
> georef2geojson VGBL42404651
> georeftats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### MGRS
``` bash
> mgrs2geojson 34TGK56063228
> mgrstats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### Tilecode
``` bash
> tilecode2geojson z23x6680749y3941729
> tilecodestats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### Maidenhead
``` bash
> maidenhead2geojson OK30is46 
> maidenheadstats # Number of cells, Average Edge Leng, Avagrae Cell Area at each resolution
```

### GARS
``` bash
> gars2geojson 574JK1918
> garsstats 574JK1918
```

### Command line for creating geocoding grids in shapefile format
``` bash
> h3grid -r 1 -o h3_1.shp (r = [0..15])
> s2grid -r 1 -o s2_1.shp (r = [0..30])
> rhealpixgrid -r 1 -o rhealpix_1.shp (r = [1..12])
> geohashgrid -r 1 -o geohash_1.shp (r = [1..12])
> gzd -o gzd.shp (Create Grid Zone Designators - used by MGRS)
> mgrsgrid -o mgrs_32648.shp -cellsize 100000 -epsg 32648 (Create MGRS Grid with cell size 100km x 100km at UTM zone 48N)  
> maidenheadgrid -r 1 -o maidenhead_1.shp (r = [1, 2, 3, 4])
```

## References

### rHEALPix
https://github.com/manaakiwhenua/rhealpixdggs-py

### DGGRID
https://github.com/sahrk/DGGRID
https://github.com/allixender/dggrid4py

### OpenEaggr DGGS
https://github.com/riskaware-ltd/open-eaggr?tab=readme-ov-file
https://github.com/riskaware-ltd/open-eaggr/blob/master/EAGGRPython/eaggr/eaggr.py
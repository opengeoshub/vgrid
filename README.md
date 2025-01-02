# Vgrid - DGGS and Cell-based Geocoding Utilities

## Installation: 
- Using pip:   
    ``` bash 
    pip install vgrid --upgrade
    ```
    
- Visit Vgrid on [PyPI](https://pypi.org/project/vgrid/)

## Demo Page:  [Vgrid Home](https://vgrid.vn)

## Usage - Vgrid CLI:
### H3
``` bash
> h32geojson 8d65b56628e46bf 
> geojson2h3 -r 11 -geojson polygon.geojson # geojson2h3 -r <res>[0..15] -geojson <point/ linestring/ polygon GeoJSON file>
> latlon2h3 10.775275567242561 106.70679737574993 13 # latlon2h3 <lat> <lon> <res> [0..15] 
> h3grid -r 11 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # h3grid -r <res> [0..15] -b <min_lon> <min_lat> <max_lon> <max_lat>
> h3stats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### S2
``` bash
> s22geojson 31752f45cc94 
> geojson2s2 -r 18 -geojson polygon.geojson # geojson2s2 -r <res>[0..30] -geojson <point/ linestring/ polygon GeoJSON file>
> latlon2s2 10.775275567242561 106.70679737574993 21 # latlon2h3 <lat> <lon> <res> [0..30]
> s2grid -r 18 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # s2grid -r <res> [0..30] -b <min_lon> <min_lat> <max_lon> <max_lat>
> s2stats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### Rhealpix
``` bash
> rhealpix2geojson R31260335553825
> geojson2rhealpix -r 11 -geojson polygon.geojson # geojson2rhealpix -r <res>[1..15] -geojson <point/ linestring/ polygon GeoJSON file>
> latlon2rhealpix 10.775275567242561 106.70679737574993 14 # latlon2rhealpix <lat> <lon> <res> [1..15]
> rhealpixgrid -r 11 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # rhealpix2grid -r <res> [0..30] -b <min_lon> <min_lat> <max_lon> <max_lat>
> rhealpixstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### OpenEAGGR ISEA4T
``` bash
> isea4t2geojson 13102313331320133331133
> geojson2isea4t -r 17 -geojson polygon.geojson # geojson2isea4t -r <res>[0..22] -geojson <point/ linestring/ polygon GeoJSON file>
> latlon2isea4t 10.775275567242561 106.70679737574993 21 # latlon2isea4t <lat> <lon> <res> [0..38]
> isea4tgrid -r 17 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea4tgrid -r <res> [0..38] -b <min_lon> <min_lat> <max_lon> <max_lat>
> isea4tstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### OpenEAGGR ISEA3H
``` bash
> latlon2isea3h 10.775275567242561 106.70679737574993 12 # latlon2isea4t <lat> <lon> <res> [0..18]
> isea3hgrid -r 0 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea3hgrid -r <res> [0..14] -b <min_lon> <min_lat> <max_lon> <max_lat>
> isea3hstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```


### OLC
``` bash
> olc2geojson 7P28QPG4+4P7
> latlon2olc 10.775275567242561 106.70679737574993 11 # latlon2olc <lat> <lon> <res> [10..15]
> olcstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### Geohash
``` bash
> geohash2geojson w3gvk1td8
> latlon2geohash 10.775275567242561 106.70679737574993 9 # latlon2geohash <lat> <lon> <res>[1..30]
> geohashgrid -r 8 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea4tgrid -r <res> [1..12] -b <min_lon> <min_lat> <max_lon> <max_lat> 1
> geohashstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### GEOREF
``` bash
> georef2geojson VGBL42404651
> latlon2georef 10.775275567242561 106.70679737574993 4 # latlon2georef <lat> <lon> <res> [0..10]
> georeftats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### MGRS
``` bash
> mgrs2geojson 34TGK56063228
> latlon2mgrs 10.775275567242561 106.70679737574993 4 # latlon2mgrs <lat> <lon> <res> [0..5]
> gzd # Create Grid Zone Designators - used by MGRS
> mgrstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### Tilecode
``` bash
> tilecode2geojson z23x6680749y3941729
> latlon2tilecode 10.775275567242561 106.70679737574993 23 # latlon2tilecode <lat> <lon> <res> [0..26]
> tilecodestats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### Maidenhead
``` bash
> maidenhead2geojson OK30is46 
> latlon2maidenhead 10.775275567242561 106.70679737574993 4 # latlon2maidenhead <lat> <lon> <res> [1..4]
> maidenheadgrid -r 4 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # maidenheadgrid -r <res> [1..4] -b <min_lon> <min_lat> <max_lon> <max_lat>
> maidenheadstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### GARS
``` bash
> gars2geojson 574JK1918
> latlon2gars 10.775275567242561 106.70679737574993 1 # latlon2gars <lat> <lon> <res> [30,15,5,1] minutes
> garsgrid -r 1 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # garsgrid -r <res> = [30,15,5,1] minutes -b <min_lon> <min_lat> <max_lon> <max_lat>
> garsstats # Number of cells, Average Edge Leng, Average Cell Area at each resolution
```

### Command line for creating DGGS and cell-based geocoding grids in GeoJSON format
``` bash
> mgrsgrid -cellsize 100000 -epsg 32648 (Create MGRS Grid with cell size 100km x 100km at UTM zone 48N)  
```

## Usage - Python code:
### Import vgrid, initialize latitude and longitude for testing:
``` python
from vgrid.utils import s2, olc, geohash, georef, mgrs, tilecode, maidenhead, gars
import h3, json
from vgrid.conversion.cell2geojson import *
from vgrid.conversion.latlon2cell import *

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
lat_lng = s2.LatLng.from_degrees(latitude, longitude)
cell_id = s2.CellId.from_lat_lng(lat_lng)
cell_id = cell_id.parent(s2_resolution)
cell_id_token= s2.CellId.to_token(cell_id)
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
gars_grid = gars.garsgrid.GARSGrid.from_latlon(latitude, longitude, gars_resolution)
gars_code = gars_grid.gars_id
print(gars_code)

data = gars2geojson(gars_code)
output_file = f'gars_{gars_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  
print(f'GeoJSON written to {output_file}')
```
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
> latlon2h3 10.775275567242561 106.70679737574993 13 # latlon2h3 <lat> <lon> <res> [0..15] 
> h32geojson 8d65b56628e46bf 
> geojson2h3 -r 11 -geojson polygon.geojson # geojson2h3 -r <res>[0..15] -geojson <point/ linestring/ polygon GeoJSON file>
> h3grid -r 11 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # h3grid -r <res> [0..15] -b <min_lon> <min_lat> <max_lon> <max_lat>
> h3stats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### S2
``` bash
> latlon2s2 10.775275567242561 106.70679737574993 21 # latlon2h3 <lat> <lon> <res> [0..30]
> s22geojson 31752f45cc94 
> geojson2s2 -r 18 -geojson polygon.geojson # geojson2s2 -r <res>[0..30] -geojson <point/ linestring/ polygon GeoJSON file>
> s2grid -r 18 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # s2grid -r <res> [0..30] -b <min_lon> <min_lat> <max_lon> <max_lat>
> s2stats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### Rhealpix
``` bash
> latlon2rhealpix 10.775275567242561 106.70679737574993 14 # latlon2rhealpix <lat> <lon> <res> [1..15]
> rhealpix2geojson R31260335553825
> geojson2rhealpix -r 11 -geojson polygon.geojson # geojson2rhealpix -r <res>[1..15] -geojson <point/ linestring/ polygon GeoJSON file>
> rhealpixgrid -r 11 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # rhealpix2grid -r <res> [0..30] -b <min_lon> <min_lat> <max_lon> <max_lat>
> rhealpixstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### OpenEAGGR ISEA4T (Windows only)
``` bash
> latlon2isea4t 10.775275567242561 106.70679737574993 21 # latlon2isea4t <lat> <lon> <res> [0..39]
> isea4t2geojson 13102313331320133331133
> geojson2isea4t -r 17 -geojson polygon.geojson # geojson2isea4t -r <res>[0..22] -geojson <point/ linestring/ polygon GeoJSON file>
> isea4tgrid -r 17 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea4tgrid -r <res> [0..25] -b <min_lon> <min_lat> <max_lon> <max_lat>
> isea4tstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### OpenEAGGR ISEA3H (Windows only)
``` bash
> latlon2isea3h 10.775275567242561 106.70679737574993 27 # latlon2isea3h <lat> <lon> <res> [0..40]
> isea3h2geojson 1327916769,-55086 
> isea3hgrid -r 20 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea3hgrid -r <res> [0..32] -b <min_lon> <min_lat> <max_lon> <max_lat>
> isea3hstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### DGGRID (Linux only)
``` bash
> latlon2dggrid 10.775275567242561 106.70679737574993 FULLER4D 10 # latlon2dggrid <lat> <lon> <dggs_type> <res> 
> dggrid2geojson 8478420 FULLER4D 10 # dggrid2geojson <lat> <lon> <dggs_type> <res> 
> geojson2dggrid -t ISEA4H -r 17 -geojson polyline.geojson # geojson2dggrid -t <DGGS Type> -r <resolution> -a <address_type> -geojson <GeoJSON path>, supporting Point/MultiPoint, LineString/MultiLineString, Polygon/MultiPolygon
> dggridgen -t ISEA3H -r 2 -a ZORDER # dggrid -t <DGGS Type> -r <res> -a<address_type>. 
> dggridstats -t FULLER3H -r 8 #dggrid -t <DGGS Type> -r <res>. 
# <DGGS Type> chosen from [SUPERFUND,PLANETRISK,ISEA3H,ISEA4H,ISEA4T,ISEA4D,ISEA43H,ISEA7H,IGEO7,FULLER3H,FULLER4H,FULLER4T,FULLER4D,FULLER43H,FULLER7H]
# <Address Type> chosen from [Q2DI,SEQNUM,INTERLEAVE,PLANE,Q2DD,PROJTRI,VERTEX2DD,AIGEN,Z3,Z3_STRING,Z7,Z7_STRING,ZORDER,ZORDER_STRING]
```

### EASE-DGGS
``` bash
> latlon2ease 10.775275567242561 106.70679737574993 6 # latlon2easedggs <lat> <lon> <res> [0..6]
> ease2geojson L6.165767.02.02.22.45.63.05
> easegrid -r 4 -b 106.6990073571 10.7628112647 106.71767427 10.778649620 # easegrid -r <res> [0..6] -b <min_lon> <min_lat> <max_lon> <max_lat>
> easestats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### POLYHEDRA GENERATOR
``` bash
> tetrahedron  # Generate Global Tetrahedron
> cube         # Generate Global Cube
> octahedron   # Generate Global Octahedron  
> icosahedron   # Generate Global Icosahedron  
``` 

### OLC
``` bash
> latlon2olc 10.775275567242561 106.70679737574993 11 # latlon2olc <lat> <lon> <res> [2,4,6,8,10,11,12,13,14,15]
> olc2geojson 7P28QPG4+4P7
> olcgrid -r 8 -b 106.6990073571 10.7628112647 106.71767427 10.778649620 # olcgrid -r <res> [2,4,6,8,10,11,12,13,14,15] -b <min_lon> <min_lat> <max_lon> <max_lat>
> olcstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### Geohash
``` bash
> latlon2geohash 10.775275567242561 106.70679737574993 9 # latlon2geohash <lat> <lon> <res>[1..10]
> geohash2geojson w3gvk1td8
> geohashgrid -r 6 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # geohashgrid -r <res> [1..10] -b <min_lon> <min_lat> <max_lon> <max_lat> 1
> geohashstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### GEOREF
``` bash
> latlon2georef 10.775275567242561 106.70679737574993 4 # latlon2georef <lat> <lon> <res> [0..5]
> georef2geojson VGBL42404651
> geohashgrid -r 2 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # geohashgrid -r <res> [0..5] -b <min_lon> <min_lat> <max_lon> <max_lat> 
> georeftats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### MGRS
``` bash
> latlon2mgrs 10.775275567242561 106.70679737574993 4 # latlon2mgrs <lat> <lon> <res> [0..5]
> mgrs2geojson 34TGK56063228
> gzd # Create Grid Zone Designators - used by MGRS
> mgrstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### Tilecode
``` bash
> latlon2tilecode 10.775275567242561 106.70679737574993 23 # latlon2tilecode <lat> <lon> <res> [0..29]
> tilecode2geojson z23x6680749y3941729
> tilecodegrid -r 20 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # tilegrid -r <res> [0..26] 
> tilecodestats # Number of cells, Cell Width, Cell Height, Cell Area at each resolution
```

### Quadkey
``` bash
> latlon2quadkey 10.775275567242561 106.70679737574993 23 # latlon2tilecode <lat> <lon> <res> [0..29]
> quadkey2geojson 13223011131020212310000
> tilegrid -r 20 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # tilegrid -r <res> [0..29] 
> tilestats # Number of cells, Cell Width, Cell Height, Cell Area at each resolution
```

### Maidenhead
``` bash
> latlon2maidenhead 10.775275567242561 106.70679737574993 4 # latlon2maidenhead <lat> <lon> <res> [1..4]
> maidenhead2geojson OK30is46 
> maidenheadgrid -r 4 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # maidenheadgrid -r <res> [1..4] -b <min_lon> <min_lat> <max_lon> <max_lat>
> maidenheadstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### GARS
``` bash
> latlon2gars 10.775275567242561 106.70679737574993 1 # latlon2gars <lat> <lon> <res> [30,15,5,1] minutes
> gars2geojson 574JK1918
> garsgrid -r 1 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # garsgrid -r <res> = [30,15,5,1] minutes -b <min_lon> <min_lat> <max_lon> <max_lat>
> garsstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

## Usage - Python code:
### Import vgrid, initialize latitude and longitude for testing:
``` python
from vgrid.utils import s2, olc, geohash, georef, mgrs, tile, maidenhead, gars
import h3, json
from vgrid.conversion.dggs2geojson import *
from vgrid.conversion.latlon2dggs import *

latitude, longitude = 10.775275567242561, 106.70679737574993
print(f'Latitude, Longitude: ({latitude}, {longitude})')
```

### H3
``` python
h3_resolution = 13 #[0..15]
h3_code = h3.latlng_to_cell(latitude, longitude, h3_resolution)
h3_decode = h3.cell_to_latlng(h3_code)
print(f'H3 code at resolution {h3_resolution}: {h3_code}')
print(f'Decode {h3_code} to WGS84: {h3_decode}')
print(f'{h3_code} to GeoJSON:\n', h32geojson(h3_code))
```

### S2
``` python
s2_resolution = 21 #[0..30]
lat_lng = s2.LatLng.from_degrees(latitude, longitude)
cell_id = s2.CellId.from_lat_lng(lat_lng)
cell_id = cell_id.parent(s2_resolution)
cell_id_token= s2.CellId.to_token(cell_id)
print(f'S2 Cell Token at resolution {s2_resolution}: {cell_id_token}')
lat_lng = cell_id.to_lat_lng() 
print(f'Decode {cell_id_token} to WGS84: {lat_lng}')
print(f'{cell_id_token} to GeoJSON:\n', s22geojson(cell_id_token))
```

### OLC
``` python
olc_resolution = 11 #[10..15]
olc_code = olc.encode(latitude, longitude, olc_resolution)
olc_decode = olc.decode(olc_code)
print(f'OLC at resolution {olc_resolution}: {olc_code}')
print(f'Decode {olc_code} to center and cell in WGS84: {olc_decode}')
print(f'{olc_code} to GeoJSON:\n', olc2geojson(olc_code))
```

### Geohash
``` python
geohash_resolution = 9 # [1..30]
geohash_code = geohash.encode(latitude, longitude, geohash_resolution)
geohash_decode = geohash.decode(geohash_code, True)
print(f'Geohash Code at resolution {geohash_resolution}: {geohash_code}')
print(f'Decode {geohash_code} to WGS84: {geohash_decode}')
print(f'{geohash_code} to GeoJSON:\n', geohash2geojson(geohash_code))
```

### GEOREF
``` python
georef_resolution = 4 # [0..10]
georef_code = georef.encode(latitude, longitude, georef_resolution)
georef_decode = georef.decode(georef_code, True)
print(f'GEOREF Code at resolution {georef_resolution}: {georef_code}')
print(f'Decode {georef_code} to WGS84: {georef_decode}')
print(f'{georef_code} to GeoJSON:\n', georef2geojson(georef_code))
```

### MGRS
``` python
mgrs_resolution = 4 # [0..5]
mgrs_code = mgrs.toMgrs(latitude, longitude, mgrs_resolution)
mgrs_code_to_wgs = mgrs.toWgs(mgrs_code)
print(f'MGRS Code at resolution {mgrs_resolution}: {mgrs_code}')
print(f'Convert {mgrs_code} to WGS84: {mgrs_code_to_wgs}')
print(f'{mgrs_code} to GeoJSON:\n', mgrs2geojson(mgrs_code))
```

### Tilecode
``` python
tile_esolution = 23  # [0..29]
tile_code = tile.latlon2tilecode(latitude, longitude, tile_esolution)
tile_encode = tile.tilecode2latlon(tile_code)
print(f'Tilecode at zoom level {tile_esolution}: {tile_code}')
print(f'Convert {tile_code} to WGS84: {tile_encode}')
print(f'{tile_code} to GeoJSON:\n', tilecode2geojson(tile_code))
```

### Maidenhead
``` python
maidenhead_resolution = 4 #[1..4]
maidenhead_code = maidenhead.toMaiden(latitude, longitude, maidenhead_resolution)
maidenGrid = maidenhead.maidenGrid(maidenhead_code)
print(f'Maidenhead Code at resolution {maidenhead_resolution}: {maidenhead_code}')
print(f'Convert {maidenhead_code} to center and cell in WGS84 = {maidenGrid}')
print(f'{maidenhead_code} to GeoJSON:\n', maidenhead2geojson(maidenhead_code))
```

### GARS
``` python
gars_resolution = 1 # [1, 5, 15, 30 minutes]
gars_grid = gars.garsgrid.GARSGrid.from_latlon(latitude, longitude, gars_resolution)
gars_code = gars_grid.gars_id
print(f'GARS code at resolution {gars_resolution}: {gars_code}')
print(f'{gars_code} to GeoJSON:\n', gars2geojson(gars_code))
```
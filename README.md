<!-- PROJECT LOGO -->
<p align="center">
    <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/vgrid192.svg" alt="Logo">
  <h2 align="center">Vgrid</h2>
  <p align="center">
    <b><i>DGGS and Cell-based Geocoding Utilities</i><b>
    <br />
  </p>
</p>
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggs.png">
</div>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Vgrid User Guide</summary>  
  <ol>
    <li>
      <a href="#vgrid-installation">Vgrid Installation</a>     
    </li>
    <li>
      <a href="#dggs-conversion">DGGS Conversion</a>
      <ul>
        <li><a href="#lat-lon-to-dggs">Lat lon to DGGS</a></li>
        <li><a href="#dggs-to-geojson">DGGS to GeoJSON</a></li>
        <li><a href="#vector-to-dggs">Vector to DGGS</a></li>
        <li><a href="#raster-to-dggs">Raster to DGGS</a></li>
        <li><a href="#dggs-compact">DGGS Compact</a></li>
        <li><a href="#dggs-expand">DGGS Expand</a></li>
      </ul>
        <li><a href="#dggs-binning">DGGS Binning</a></li>
        <li><a href="#dggs-resampling">DGGS Resampling</a></li>
        <li><a href="#dggs-generator">DGGS Generator</a></li>
    </li>   
  </ol>
</details>

## Vgrid Installation
- Using pip:   
    ``` bash 
    pip install vgrid --upgrade
    ```
- Visit Vgrid on [PyPI](https://pypi.org/project/vgrid/)

- Demo Page:  [Vgrid Home](https://vgrid.vn)

## ***The following Vgrid user guide is intended for use in a Command Line Interface (CLI) environment.***

## DGGS Conversion

### Lat lon to DGGS

Convert lat, long in WGS84 CRS to DGGS cellID (H3, S2, rHEALPix, OpenEaggr ISEA4T and ISEA3H (Windows only), EASE-DGGS, QTM, OLC/ OpenLocationCode/ Google Plus Code, Geohash, GEOREF, MGRS, Tilecode, Quadkey, Maidenhead, GARS).
``` bash
> latlon2h3 10.775276 106.706797 13 # latlon2h3 <lat> <lon> <resolution> [0..15] 
> latlon2s2 10.775276 106.706797 21 # latlon2s2 <lat> <lon> <resolution> [0..30]
> latlon2rhealpix 10.775276 106.706797 14 # latlon2rhealpix <lat> <lon> <resolution> [1..15]
> latlon2isea4t 10.775276 106.706797 21 # latlon2isea4t <lat> <lon> <resolution> [0..39]
> latlon2isea3h 10.775276 106.706797 27 # latlon2isea3h <lat> <lon> <resolution> [0..40]
> latlon2ease 10.775276 106.706797 6 # latlon2easedggs <lat> <lon> <resolution> [0..6]
> latlon2dggrid 10.775276 106.706797 FULLER4D 10 # Linux only: latlon2dggrid <lat> <lon> <dggs_type> <resolution> 
> latlon2qtm 10.775276 106.706797 11  # latlon2qtm <lat> <lon> <resolution> [1..24]
> latlon2olc 10.775276 106.706797 11 # latlon2olc <lat> <lon> <resolution> [2,4,6,8,10,11..15]
> latlon2geohash 10.775276 106.706797 9 # latlon2geohash <lat> <lon> <resolution>[1..10]
> latlon2georef 10.775276 106.706797 4 # latlon2georef <lat> <lon> <resolution> [0..5]
> latlon2mgrs 10.775276 106.706797 4 # latlon2mgrs <lat> <lon> <resolution> [0..5]
> latlon2tilecode 10.775276 106.706797 23 # latlon2tilecode <lat> <lon> <resolution> [0..29]
> latlon2quadkey 10.775276 106.706797 23 # latlon2quadkey <lat> <lon> <resolution> [0..29]
> latlon2maidenhead 10.775276 106.706797 4 # latlon2maidenhead <lat> <lon> <resolution> [1..4]
> latlon2gars 10.775276 106.706797 1 # latlon2gars <lat> <lon> <resolution> [30,15,5,1] minutes
```

### DGGS to GeoJSON

Convert DGGS cell ID to GeoJSON.
``` bash
> h32geojson 8d65b56628e46bf 
> s22geojson 31752f45cc94 
> rhealpix2geojson R31260335553825
> isea4t2geojson 13102313331320133331133 # Windows only
> isea3h2geojson 1327916769,-55086 # Windows only
> ease2geojson L6.165767.02.02.22.45.63.05
> dggrid2geojson 8478420 FULLER4D 10 # Linux only:  dggrid2geojson <DGGRID code in SEQNUM format> <DGGS Type> <resolution>
dggs_types = (
    'CUSTOM',  # parameters will be specified manually
    'SUPERFUND', # Superfund_500m grid
    'PLANETRISK',
    'ISEA3H', # ISEA projection with hexagon cells and an aperture of 3
    'ISEA4H', # ISEA projection with hexagon cells and an aperture of 4
    'ISEA4T',  # ISEA projection with triangle cells and an aperture of 4
    'ISEA4D', # ISEA projection with diamond cells and an aperture of 4
    'ISEA43H', # ISEA projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
    'ISEA7H', # ISEA projection with hexagon cells and an aperture of 7
    'IGEO7', # ISEA projection with hexagon cells and an aperture of 7 and Z7 address type
    'FULLER3H', # FULLER projection with hexagon cells and an aperture of 3
    'FULLER4H', # FULLER projection with hexagon cells and an aperture of 4
    'FULLER4T', # FULLER projection with triangle cells and an aperture of 4
    'FULLER4D', # FULLER projection with diamond cells and an aperture of 4
    'FULLER43H', # FULLER projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
    'FULLER7H', # FULLER projection with hexagon cells and an aperture of 7
)
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
```

### Vector to DGGS
Convert Vector layers (Point/ Multipoint, Linestring/ Multilinestring, polygon/ Multipolygon) in GeoJSON to DGGS.

- **Uncompact:**
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/vertor2dggs_uncompact.png">
</div>

- **Compact:**
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/vertor2dggs_compact.png">
</div>

``` bash
> geojson2h3 -r 11 -geojson polygon.geojson # geojson2h3 -r <resolution>[0..15] -geojson <GeoJSON file> -compact [optional]
> geojson2s2 -r 18 -geojson polygon.geojson -compact # geojson2s2 -r <resolution>[0..30] -geojson <GeoJSON file> -compact [optional]
> geojson2rhealpix -r 11 -geojson polygon.geojson # geojson2rhealpix -r <resolution>[1..15] -geojson <GeoJSON file> -compact [optional]
> geojson2isea4t -r 17 -geojson polygon.geojson # geojson2isea4t -r <resolution>[0..25] -geojson <GeoJSON file> -compact [optional]
> geojson2isea3h -r 18 -geojson polygon.geojson  # geojson2isea3h -r <resolution>[1..32] -geojson <GeoJSON file> -compact [optional]
> geojson2ease -r 4 -geojson polygon.geojson  # geojson2ease -r <resolution>[0..6] -geojson <GeoJSON file> -compact [optional]
> geojson2dggrid -t ISEA4H -r 17 -geojson polyline.geojson # Linux only: geojson2dggrid -t <DGGS Type> -r <resolution> -a <address_type, default is SEQNUM> -geojson <GeoJSON path>
> geojson2qtm -r 18 -geojson polygon.geojson  # geojson2qtm -r <resolution>[1..24] -geojson <GeoJSON file> -compact [optional]
> geojson2olc -r 10 -geojson polygon.geojson # geojson2olc -r <resolution>[2,4,6,8,10,11..15] -geojson <GeoJSON file> -compact [optional]
> geojson2geohash -r 7 -geojson polygon.geojson # geojson2geohash -r <resolution>[1..10] -geojson <GeoJSON file> -compact [optional]
> geojson2mgrs -r 3 -geojson polygon.geojson # geojson2mgrs -r <resolution>[0..5] -geojson <GeoJSON file>
> geojson2tilecode -r 18 -geojson polygon.geojson # geojson2tilecode -r <resolution>[0..29] -geojson <GeoJSON file> -compact [optional]
> geojson2quadkey -r 18 -geojson polygon.geojson # geojson2quadkey -r <resolution>[0..29] -geojson <GeoJSON file> -compact [optional]
```
### Raster to DGGS

Convert raster layers in geographic CRS to DGGS.
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/raster2dggs_h3.png">
</div>

``` bash
> raster2h3 -raster raster.tif # raster2h3 -raster <raster in geographic CRS> -r <resolution>[0..15] [optional, defaults to the H3 resolution nearest to the raster's cell size]
> raster2s2 -raster raster.tif # raster2s2 -raster <raster in geographic CRS> -r <resolution>[0..24] [optional, defaults to the S2 resolution nearest to the raster's cell size]
> raster2rhealpix -raster raster.tif # raster2rhealpix -raster <raster in geographic CRS> -r <resolution>[0..15] [optional, defaults to the rHEALPix resolution nearest to the raster's cell size]
> raster2isea4t -raster raster.tif # Windows only: raster2isea4t -raster <raster in geographic CRS> -r <resolution>[0..23] [optional, defaults to the ISEA4T resolution nearest to the raster's cell size]
> raster2qtm -raster raster.tif # raster2qtm -raster <raster in geographic CRS> -r <resolution>[1..24] [optional, defaults to the QTM resolution nearest to the raster's cell size]
> raster2olc -raster raster.tif # raster2olc -raster <raster in geographic CRS> -r <resolution>[10..12] [optional, defaults to the OLC resolution nearest to the raster's cell size]
> raster2geohash -raster raster.tif # raster2geohash -raster <raster in geographic CRS> -r <resolution>[1..10] [optional, defaults to the Geohash resolution nearest to the raster's cell size]
> raster2tilecode -raster raster.tif # raster2tilecode -raster <raster in geographic CRS> -r <resolution>[0..26] [optional, defaults to the Tilecode resolution nearest to the raster's cell size]
> raster2quadkey -raster raster.tif # raster2quadkey -raster <raster in geographic CRS> -r <resolution>[0..26] [optional, defaults to the Quadkey resolution nearest to the raster's cell size]
```

### DGGS Compact
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggscompact_isea4t.png">
</div>

``` bash
> h3compact -geojson h3.geojson -cellid h3  # h3expand -geojson <H3 in GeoJSON> -cellid [optional, 'h3' by default]
> s2compact -geojson s2.geojson -cellid s2  # h3expand -geojson <S2 in GeoJSON> -cellid [optional, 's2' by default]
> rhealpixcompact -geojson rhealpix.geojson -cellid rhealpix  # h3expand -geojson <rHEALPix in GeoJSON> -cellid [optional, 'rhealpix' by default]
```
### DGGS Expand
<div align="center">
  <img src="https://raw.githubusercontent.com/thangqd/vgridtools/main/images/readme/dggsexpand_isea4t.png">
</div>

``` bash
> h3expand -geojson h3_11.geojson -r 12 -cellid h3 # h3expand -geojson <H3 in smaller resolition> -r <higher resolution>[0..15] -cellid [optional, 'h3' by default]
```

## DGGS Binning
### Binning point layer to DGGS.
<div align="center">
  <img src="images/readme/dggsbinning.png">
</div>

<div align="center">
  <img src="images/readme/dggsbinning_h3.png">
</div>

## DGGS Resampling
<div align="center">
  <img src="images/readme/dggsresampling.png">
</div>

<div align="center">
  <img src="images/readme/dggsresampling_h32s2.png">
</div>

## DGGS Generator
<div align="center">
  <img src="images/readme/dggsgenerator.png">
</div>

<div align="center">
  <img src="images/readme/dggsgenerator_h3.png">
</div>


### DGGS Resampling
``` bash
> resample -geojson vn_pop_h3_6.geojson -fromdggs h3 -todggs s2 -resamplefield population # resample -geojson <Input DGGS> -fromdggs <Input DGGS type> -todggs <Output DGGS Type> -resamplefield [Optional, Numeric field for resampling]
```

### H3
``` bash
> csv2h3 h3.csv  # Convert CSV with 'h3' column to H3 cells.
> h3bin -point point.geojson -r 8 -stats count -field numeric_field -category group # h3bin -point <point GeoJSON file> -r <resolutin[0..15]> -stats <count, min, max, sum, mean, median, std, var, range, minority, majority, variety> -field [Optional, numeric field to compute statistics] -category [optional, category field for grouping] 
> h3grid -r 11 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # h3grid -r <resolution> [0..15] -b <min_lon> <min_lat> <max_lon> <max_lat>
> h3stats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### S2
``` bash
> s2grid -r 18 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # s2grid -r <resolution> [0..30] -b <min_lon> <min_lat> <max_lon> <max_lat>
> s2stats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### Rhealpix
``` bash
> rhealpixgrid -r 11 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # rhealpix2grid -r <resolution> [0..30] -b <min_lon> <min_lat> <max_lon> <max_lat>
> rhealpixstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### OpenEAGGR ISEA4T (Windows only)
``` bash
> isea4tgrid -r 17 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea4tgrid -r <resolution> [0..25] -b <min_lon> <min_lat> <max_lon> <max_lat>
> isea4tstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### OpenEAGGR ISEA3H (Windows only)
``` bash
> isea3hgrid -r 20 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # isea3hgrid -r <resolution> [0..32] -b <min_lon> <min_lat> <max_lon> <max_lat>
> isea3hstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### DGGRID (Linux only)
``` bash

> dggridgen -t ISEA3H -r 2 -a ZORDER # dggrid -t <DGGS Type> -r <resolution> -a<address_type>. 
> dggridstats -t FULLER3H -r 8 #dggrid -t <DGGS Type> -r <resolution>. 
# <DGGS Type> chosen from [SUPERFUND,PLANETRISK,ISEA3H,ISEA4H,ISEA4T,ISEA4D,ISEA43H,ISEA7H,IGEO7,FULLER3H,FULLER4H,FULLER4T,FULLER4D,FULLER43H,FULLER7H]
# <Address Type> chosen from [Q2DI,SEQNUM,INTERLEAVE,PLANE,Q2DD,PROJTRI,VERTEX2DD,AIGEN,Z3,Z3_STRING,Z7,Z7_STRING,ZORDER,ZORDER_STRING]
```

### EASE-DGGS
``` bash
> easegrid -r 4 -b 106.6990073571 10.7628112647 106.71767427 10.778649620 # easegrid -r <resolution> [0..6] -b <min_lon> <min_lat> <max_lon> <max_lat>
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
> olcgrid -r 8 -b 106.6990073571 10.7628112647 106.71767427 10.778649620 # olcgrid -r <resolution> [2,4,6,8,10,11,12,13,14,15] -b <min_lon> <min_lat> <max_lon> <max_lat>
> olcstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### Geohash
``` bash
> geohashgrid -r 6 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # geohashgrid -r <resolution> [1..10] -b <min_lon> <min_lat> <max_lon> <max_lat> 1
> geohashstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### GEOREF
``` bash
> georef2geojson VGBL42404651
> geohashgrid -r 2 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # geohashgrid -r <resolution> [0..5] -b <min_lon> <min_lat> <max_lon> <max_lat> 
> georeftats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### MGRS
``` bash
> mgrstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### Tilecode
``` bash
> tilecodegrid -r 20 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # tilegrid -r <resolution> [0..26] 
> tilecodestats # Number of cells, Cell Width, Cell Height, Cell Area at each resolution
```

### Quadkey
``` bash
> tilegrid -r 20 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # tilegrid -r <resolution> [0..29] 
> tilestats # Number of cells, Cell Width, Cell Height, Cell Area at each resolution
```

### Maidenhead
``` bash
> maidenheadgrid -r 4 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # maidenheadgrid -r <resolution> [1..4] -b <min_lon> <min_lat> <max_lon> <max_lat>
> maidenheadstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```

### GARS
``` bash
> garsgrid -r 1 -b 106.6990073571 10.7628112647 106.71767427 10.7786496202 # garsgrid -r <resolution> = [30,15,5,1] minutes -b <min_lon> <min_lat> <max_lon> <max_lat>
> garsstats # Number of cells, Avg Edge Length, Avg Cell Area at each resolution
```
# python setup.py sdist bdist_wheel
# twine upload dist/*
import os
import shutil
from setuptools import setup, find_packages

requirements = [
    'tqdm~=4.66.2',
    'shapely~=2.0.1',
    'protobuf~=5.26.1',
    'fiona',
    'pyproj',
    'pyclipper~=1.3.0',
    'h3~=4.1.1',
    'geopandas',
    'scipy',
    'future',
    'texttable',
    'rasterio'
    ],

def clean_build():
    build_dir = 'build'
    dist_dir = 'dist'
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
    if os.path.exists(dist_dir):
        shutil.rmtree(dist_dir)

clean_build()

setup(
    name='vgrid',
    version='1.3.5',
    author = 'Thang Quach',
    author_email= 'quachdongthang@gmail.com',
    url='https://github.com/thangqd/vgrid',
    description='Vgrid - DGGS and Cell-based Geocoding Utilites',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    requires_python=">=3.0",
    packages=find_packages(),
    include_package_data=True,  # Include package data specified in MANIFEST.in
    entry_points={
        'console_scripts': [  
            # Latlon to DGGS Cell
            'latlon2h3 = vgrid.conversion.latlon2dggs:latlon2h3_cli',  
            'latlon2s2 = vgrid.conversion.latlon2dggs:latlon2s2_cli',  
            'latlon2rhealpix = vgrid.conversion.latlon2dggs:latlon2rhealpix_cli',  
            'latlon2isea4t = vgrid.conversion.latlon2dggs:latlon2isea4t_cli',  
            'latlon2isea3h = vgrid.conversion.latlon2dggs:latlon2isea3h_cli',  
            'latlon2ease = vgrid.conversion.latlon2dggs:latlon2ease_cli',
            
            'latlon2dggrid = vgrid.conversion.latlon2dggs:latlon2dggrid_cli',
            
            'latlon2qtm = vgrid.conversion.latlon2dggs:latlon2qtm_cli',

            'latlon2olc = vgrid.conversion.latlon2dggs:latlon2olc_cli',  
            'latlon2geohash = vgrid.conversion.latlon2dggs:latlon2geohash_cli',  
            'latlon2georef = vgrid.conversion.latlon2dggs:latlon2georef_cli',  
            'latlon2mgrs = vgrid.conversion.latlon2dggs:latlon2mgrs_cli',  
            'latlon2tilecode = vgrid.conversion.latlon2dggs:latlon2tilecode_cli',  
            'latlon2quadkey = vgrid.conversion.latlon2dggs:latlon2quadkey_cli',             
            'latlon2maidenhead = vgrid.conversion.latlon2dggs:latlon2maidenhead_cli',  
            'latlon2gars = vgrid.conversion.latlon2dggs:latlon2gars_cli',  

            # DGGS ID to GeoJSON
            'h32geojson = vgrid.conversion.dggs2geojson:h32geojson_cli',  
            's22geojson = vgrid.conversion.dggs2geojson:s22geojson_cli',  
            'rhealpix2geojson = vgrid.conversion.dggs2geojson:rhealpix2geojson_cli',  
            'isea4t2geojson = vgrid.conversion.dggs2geojson:isea4t2geojson_cli',  
            
            'isea3h2geojson = vgrid.conversion.dggs2geojson:isea3h2geojson_cli',  

            'dggrid2geojson = vgrid.conversion.dggs2geojson:dggrid2geojson_cli',  

            'ease2geojson = vgrid.conversion.dggs2geojson:ease2geojson_cli',  
            'qtm2geojson = vgrid.conversion.dggs2geojson:qtm2geojson_cli',  

            'olc2geojson = vgrid.conversion.dggs2geojson:olc2geojson_cli',
            'geohash2geojson = vgrid.conversion.dggs2geojson:geohash2geojson_cli',  
            'georef2geojson = vgrid.conversion.dggs2geojson:georef2geojson_cli',  
            'mgrs2geojson = vgrid.conversion.dggs2geojson:mgrs2geojson_cli', 
            'tilecode2geojson = vgrid.conversion.dggs2geojson:tilecode2geojson_cli',  
            'quadkey2geojson = vgrid.conversion.dggs2geojson:quadkey2geojson_cli',  

            'maidenhead2geojson = vgrid.conversion.dggs2geojson:maidenhead2geojson_cli',  
            'gars2geojson = vgrid.conversion.dggs2geojson:gars2geojson_cli',  
            
            # GeoJSON to DGGS
            'geojson2h3 = vgrid.conversion.geojson2dggs.geojson2h3:main',
            'geojson2s2 = vgrid.conversion.geojson2dggs.geojson2s2:main',
            'geojson2rhealpix = vgrid.conversion.geojson2dggs.geojson2rhealpix:main',
            'geojson2isea4t = vgrid.conversion.geojson2dggs.geojson2isea4t:main',
            'geojson2isea3h = vgrid.conversion.geojson2dggs.geojson2isea3h:main',
            'geojson2ease = vgrid.conversion.geojson2dggs.geojson2ease:main',

            'geojson2dggrid = vgrid.conversion.geojson2dggs.geojson2dggrid:main',
            'geojson2qtm = vgrid.conversion.geojson2dggs.geojson2qtm:main',
            
            'geojson2olc = vgrid.conversion.geojson2dggs.geojson2olc:main',
            'geojson2geohash = vgrid.conversion.geojson2dggs.geojson2geohash:main',
            # 'geojson2georef = vgrid.conversion.geojson2dggs.geojson2georef:main',
            'geojson2mgrs = vgrid.conversion.geojson2dggs.geojson2mgrs:main',
            'geojson2tilecode = vgrid.conversion.geojson2dggs.geojson2tilecode:main',
            'geojson2quadkey = vgrid.conversion.geojson2dggs.geojson2quadkey:main',

            # DGGS compact/ expand
            'h3compact =  vgrid.conversion.dggscompact:h3compact_cli',
            'h3expand =  vgrid.conversion.dggscompact:h3expand_cli',
            's2compact =  vgrid.conversion.dggscompact:s2compact_cli',
            's2expand =  vgrid.conversion.dggscompact:s2expand_cli',
            'rhealpixcompact = vgrid.conversion.dggscompact:rhealpixcompact_cli',
            'rhealpixexpand = vgrid.conversion.dggscompact:rhealpixexpand_cli',
            'isea4tcompact = vgrid.conversion.dggscompact:isea4tcompact_cli',
            'isea4texpand = vgrid.conversion.dggscompact:isea4texpand_cli',
            'isea3hcompact = vgrid.conversion.dggscompact:isea3hcompact_cli',
            'isea3hexpand = vgrid.conversion.dggscompact:isea3hexpand_cli',
            'easecompact = vgrid.conversion.dggscompact:easecompact_cli',
            'easeexpand = vgrid.conversion.dggscompact:easeexpand_cli',
            'qtmcompact = vgrid.conversion.dggscompact:qtmcompact_cli',
            'qtmexpand = vgrid.conversion.dggscompact:qtmexpand_cli',
            

            'olccompact = vgrid.conversion.dggscompact:olccompact_cli',
            'olcexpand = vgrid.conversion.dggscompact:olcexpand_cli',
            'geohashcompact = vgrid.conversion.dggscompact:geohashcompact_cli',
            'geohashexpand = vgrid.conversion.dggscompact:geohashexpand_cli',
            'tilecodecompact = vgrid.conversion.dggscompact:tilecodecompact_cli',
            'tilecodeexpand = vgrid.conversion.dggscompact:tilecodeexpand_cli',
            'quadkeycompact = vgrid.conversion.dggscompact:quadkeycompact_cli',
            'quadkeyexpand = vgrid.conversion.dggscompact:quadkeyexpand_cli',


            # Raster to DGGS
            'raster2h3 = vgrid.conversion.raster2dggs.raster2h3:main',
            'raster2s2 = vgrid.conversion.raster2dggs.raster2s2:main',
            'raster2rhealpix = vgrid.conversion.raster2dggs.raster2rhealpix:main',
            'raster2isea4t = vgrid.conversion.raster2dggs.raster2isea4t:main',
            'raster2qtm = vgrid.conversion.raster2dggs.raster2qtm:main',
            'raster2geohash = vgrid.conversion.raster2dggs.raster2geohash:main',
            'raster2olc = vgrid.conversion.raster2dggs.raster2olc:main',
            'raster2tilecode = vgrid.conversion.raster2dggs.raster2tilecode:main',
            'raster2quadkey = vgrid.conversion.raster2dggs.raster2quadkey:main',

            # CSV to DGGS Cell 
            'csv2h3 =  vgrid.conversion.csv2dggs:csv2h3',
            'csv2s2 =  vgrid.conversion.csv2dggs:csv2s2',
            'csv2rhealpix =  vgrid.conversion.csv2dggs:csv2rhealpix',
            'csv2isea4t =  vgrid.conversion.csv2dggs:csv2isea4t',
            'csv2isea3h =  vgrid.conversion.csv2dggs:csv2isea3h',
            'csv2ease =  vgrid.conversion.csv2dggs:csv2ease',
            'csv2qtm =  vgrid.conversion.csv2dggs:csv2qtm',

            'csv2olc =  vgrid.conversion.csv2dggs:csv2olc',
            'csv2geohash =  vgrid.conversion.csv2dggs:csv2geohash',
            'csv2georef =  vgrid.conversion.csv2dggs:csv2georef',
            'csv2mgrs =  vgrid.conversion.csv2dggs:csv2mgrs',
            'csv2tilecode =  vgrid.conversion.csv2dggs:csv2tilecode',
            'csv2quadkey =  vgrid.conversion.csv2dggs:csv2quadkey',
            'csv2maidenhead =  vgrid.conversion.csv2dggs:csv2maidenhead',
            'csv2gars =  vgrid.conversion.csv2dggs:csv2gars',
            
            # Data binning
            'h3bin =  vgrid.binning.h3bin:main',
            's2bin =  vgrid.binning.s2bin:main',
            'rhealpixbin =  vgrid.binning.rhealpixbin:main',
            'isea4tbin =  vgrid.binning.isea4tbin:main',
            'qtmbin =  vgrid.binning.qtmbin:main',
            
            'olcbin =  vgrid.binning.olcbin:main',
            'geohashbin =  vgrid.binning.geohashbin:main',
            'tilecodebin =  vgrid.binning.tilecodebin:main',
            'quadkeybin =  vgrid.binning.quadkeybin:main',

            'polygonbin =  vgrid.binning.polygonbin:main',

            # Resampling
            'dggsresample =  vgrid.resampling.dggsresample:main',
            
            # Grid Generator
            'h3grid = vgrid.generator.h3grid:main',
            's2grid = vgrid.generator.s2grid:main',
            'rhealpixgrid = vgrid.generator.rhealpixgrid:main',
            'isea4tgrid = vgrid.generator.isea4tgrid:main',
            'isea3hgrid = vgrid.generator.isea3hgrid:main',
            'easegrid = vgrid.generator.easegrid:main', # need to be checked
            'dggridgen = vgrid.generator.dggridgen:main',
            'qtmgrid = vgrid.generator.qtmgrid:main',

            'olcgrid = vgrid.generator.olcgrid:main',
            'geohashgrid = vgrid.generator.geohashgrid:main',    
            'georefgrid = vgrid.generator.georefgrid:main',           
            'gzd = vgrid.generator.gzd:main',  
            'mgrsgrid = vgrid.generator.mgrsgrid:main',
            'tilecodegrid = vgrid.generator.tilecodegrid:main', 
            'quadkeygrid = vgrid.generator.quadkeygrid:main', 
            'maidenheadgrid = vgrid.generator.maidenheadgrid:main',        
            'garsgrid = vgrid.generator.garsgrid:main',        
            
            # Polyhedra Generator  
            'tetrahedron = vgrid.generator.polyhedra.tetrahedron:main',   
            'cube = vgrid.generator.polyhedra.cube:main',        
            'octahedron = vgrid.generator.polyhedra.octahedron:main',        
            'hexagon = vgrid.generator.polyhedra.hexagon:main',  
            'fuller_icosahedron = vgrid.generator.polyhedra.fuller_icosahedron:main',        
            'rhombic_icosahedron = vgrid.generator.polyhedra.rhombic_icosahedron:main',        

             # Grid Stats
            'h3stats = vgrid.stats.h3stats:main',
            's2stats = vgrid.stats.s2stats:main',
            'rhealpixstats = vgrid.stats.rhealpixstats:main',
            'isea4tstats = vgrid.stats.isea4tstats:main',
            'isea3hstats = vgrid.stats.isea3hstats:main',
            'easestats = vgrid.stats.easestats:main',
            'qtmstats = vgrid.stats.qtmstats:main',

            'dggridstats = vgrid.stats.dggridstats:main',

            'olcstats = vgrid.stats.olcstats:main',
            'geohashstats = vgrid.stats.geohashstats:main',
            'georefstats = vgrid.stats.georefstats:main',
            'mgrsstats = vgrid.stats.mgrsstats:main',
            'tilecodestats = vgrid.stats.tilecodestats:main',
            'quadkeystats = vgrid.stats.quadkeystats:main',

            'maidenheadstats = vgrid.stats.maidenheadstats:main',
            'garsstats = vgrid.stats.garsstats:main',

            # DGGRID Corrections
            'dggridfixcontent = vgrid.correction.dggridfixcontent:main',
            'dggridfixgeom = vgrid.correction.dggridfixgeom:main',       
            'dggridfixgeom2 = vgrid.correction.dggridfixgeom2:main',            
            'dggridfix = vgrid.correction.dggridfix:main',    
        ],
    },    

    install_requires=requirements,    
    classifiers=[
        'Programming Language :: Python :: 3',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: GIS',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)

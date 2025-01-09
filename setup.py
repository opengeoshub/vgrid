# python setup.py sdist bdist_wheel
# twine upload dist/*
import os
import shutil
from setuptools import setup, find_packages

requirements = [
    'tqdm~=4.66.2',
    'shapely~=2.0.1',
    'protobuf~=5.26.1',
    'fiona~=1.10.0',
    'pyproj',
    'pyclipper~=1.3.0',
    'h3~=4.1.1',
    'geopandas',
    'scipy',
    'future',
    'texttable'
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
    version='1.2.2',
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
            # Latlon to Cell
            'latlon2h3 = vgrid.conversion.latlon2cell:latlon2h3_cli',  
            'latlon2s2 = vgrid.conversion.latlon2cell:latlon2s2_cli',  
            'latlon2rhealpix = vgrid.conversion.latlon2cell:latlon2rhealpix_cli',  
            'latlon2isea4t = vgrid.conversion.latlon2cell:latlon2isea4t_cli',  
            'latlon2isea3h = vgrid.conversion.latlon2cell:latlon2isea3h_cli',  
            'latlon2ease = vgrid.conversion.latlon2cell:latlon2ease_cli',
            
            'latlon2olc = vgrid.conversion.latlon2cell:latlon2olc_cli',  
            'latlon2geohash = vgrid.conversion.latlon2cell:latlon2geohash_cli',  
            'latlon2georef = vgrid.conversion.latlon2cell:latlon2georef_cli',  
            'latlon2mgrs = vgrid.conversion.latlon2cell:latlon2mgrs_cli',  
            'latlon2tilecode = vgrid.conversion.latlon2cell:latlon2tilecode_cli',  
            'latlon2quadkey = vgrid.conversion.latlon2cell:latlon2quadkey_cli',             
            'latlon2maidenhead = vgrid.conversion.latlon2cell:latlon2maidenhead_cli',  
            'latlon2gars = vgrid.conversion.latlon2cell:latlon2gars_cli',  

            # Cell to GeoJSON
            'h32geojson = vgrid.conversion.cell2geojson:h32geojson_cli',  
            's22geojson = vgrid.conversion.cell2geojson:s22geojson_cli',  
            'rhealpix2geojson = vgrid.conversion.cell2geojson:rhealpix2geojson_cli',  
            'isea4t2geojson = vgrid.conversion.cell2geojson:isea4t2geojson_cli',  
            
            'isea3h2geojson = vgrid.conversion.cell2geojson:isea3h2geojson_cli',  
            'isea3hparent2geojson = vgrid.utils.eaggr.conversion.parent2geojson:parent2geojson_cli',  
            'isea3hsiblings2geojson = vgrid.utils.eaggr.conversion.siblings2geojson:siblings2geojson_cli',  
            'isea3hchildren2geojson = vgrid.utils.eaggr.conversion.children2geojson:children2geojson_cli',  
            
            'ease2geojson = vgrid.conversion.cell2geojson:ease2geojson_cli',  
            
            'olc2geojson = vgrid.conversion.cell2geojson:olc2geojson_cli',
            'geohash2geojson = vgrid.conversion.cell2geojson:geohash2geojson_cli',  
            'georef2geojson = vgrid.conversion.cell2geojson:georef2geojson_cli',  
            'mgrs2geojson = vgrid.conversion.cell2geojson:mgrs2geojson_cli', 
            'tilecode2geojson = vgrid.conversion.cell2geojson:tilecode2geojson_cli',  
            'quadkey2geojson = vgrid.conversion.cell2geojson:quadkey2geojson_cli',  

            'maidenhead2geojson = vgrid.conversion.cell2geojson:maidenhead2geojson_cli',  
            'gars2geojson = vgrid.conversion.cell2geojson:gars2geojson_cli',  
            
            # GeoJSON to Cell
            'geojson2h3 = vgrid.conversion.geojson2h3:main',
            'geojson2s2 = vgrid.conversion.geojson2s2:main',
            'geojson2rhealpix = vgrid.conversion.geojson2rhealpix:main',
            'geojson2isea4t = vgrid.conversion.geojson2isea4t:main',
            'geojson2isea3h = vgrid.conversion.geojson2isea3h:main',

            # Grid Generator
            'h3grid = vgrid.generator.h3grid:main',
            's2grid = vgrid.generator.s2grid:main',
            'rhealpixgrid = vgrid.generator.rhealpixgrid:main',
            'isea4tgrid = vgrid.generator.isea4tgrid:main',
            'isea3hgrid = vgrid.generator.isea3hgrid:main',
            'easegrid = vgrid.generator.easegrid:main',

            'olcgrid = vgrid.generator.olcgrid:main',
            'geohashgrid = vgrid.generator.geohashgrid:main',    
            'georefgrid = vgrid.generator.georefgrid:main',           
            'gzd = vgrid.generator.gzd:main',  
            'mgrsgrid = vgrid.generator.mgrsgrid:main',
            'tilegrid = vgrid.generator.tilegrid:main', 
            'maidenheadgrid = vgrid.generator.maidenheadgrid:main',        
            'garsgrid = vgrid.generator.garsgrid:main',        
   
            # Grid Stats
            'h3stats = vgrid.stats.h3stats:main',
            's2stats = vgrid.stats.s2stats:main',
            'rhealpixstats = vgrid.stats.rhealpixstats:main',
            'isea4tstats = vgrid.stats.isea4tstats:main',
            'isea3hstats = vgrid.stats.isea3hstats:main',
            'easestats = vgrid.stats.easestats:main',

            'olcstats = vgrid.stats.olcstats:main',
            'geohashstats = vgrid.stats.geohashstats:main',
            'georefstats = vgrid.stats.georefstats:main',
            'mgrsstats = vgrid.stats.mgrsstats:main',
            'tilestats = vgrid.stats.tilestats:main',
            
            'maidenheadstats = vgrid.stats.maidenheadstats:main',
            'garsstats = vgrid.stats.garsstats:main',

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

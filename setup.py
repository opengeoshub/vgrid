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
    'pandas~=2.0.3',
    'geopandas',
    'scipy',
    'future',
    'rhealpixdggs'
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
    version='1.1.26',
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
            # Geocode to GeoJSON
            'h32geojson = vgrid.geocode.geocode2geojson:h32geojson_cli',  
            's22geojson = vgrid.geocode.geocode2geojson:s22geojson_cli',  
            'rhealpix2geojson = vgrid.geocode.geocode2geojson:rhealpix2geojson_cli',  
            'eaggrisea4t2geojson = vgrid.geocode.geocode2geojson:eaggrisea4t2geojson_cli',  
            'eaggrisea3h2geojson = vgrid.geocode.geocode2geojson:eaggrisea3h2geojson_cli',  

            'olc2geojson = vgrid.geocode.geocode2geojson:olc2geojson_cli',
            'geohash2geojson = vgrid.geocode.geocode2geojson:geohash2geojson_cli',  
            'georef2geojson = vgrid.geocode.geocode2geojson:georef2geojson_cli',  
            'mgrs2geojson = vgrid.geocode.geocode2geojson:mgrs2geojson_cli', 
            'tilecode2geojson = vgrid.geocode.geocode2geojson:tilecode2geojson_cli',  

            'maidenhead2geojson = vgrid.geocode.geocode2geojson:maidenhead2geojson_cli',  
            'gars2geojson = vgrid.geocode.geocode2geojson:gars2geojson_cli',  
            
            # Latlon to Code
            'latlon2h3 = vgrid.geocode.latlon2geocode:latlon2h3_cli',  
            'latlon2s2 = vgrid.geocode.latlon2geocode:latlon2s2_cli',  
            'latlon2rhealpix = vgrid.geocode.latlon2geocode:latlon2rhealpix_cli',  
            'latlon2eaggrisea4t = vgrid.geocode.latlon2geocode:latlon2eaggrisea4t_cli',  
            'latlon2olc = vgrid.geocode.latlon2geocode:latlon2olc_cli',  
            'latlon2geohash = vgrid.geocode.latlon2geocode:latlon2geohash_cli',  
            'latlon2georef = vgrid.geocode.latlon2geocode:latlon2georef_cli',  
            'latlon2mgrs = vgrid.geocode.latlon2geocode:latlon2mgrs_cli',  
            'latlon2tilecode = vgrid.geocode.latlon2geocode:latlon2tilecode_cli',  
            'latlon2maidenhead = vgrid.geocode.latlon2geocode:latlon2maidenhead_cli',  
            'latlon2gars = vgrid.geocode.latlon2geocode:latlon2gars_cli',  

            # geojson2grid
            'geojson2rhealpix = vgrid.grid.geojson2grid.geojson2rhealpix:main',
            'geojson2eaggrisea4t = vgrid.grid.geojson2grid.geojson2eaggrisea4t:main',
            'geojson2h3 = vgrid.grid.geojson2grid.geojson2h3:main',

            # Generate Geocode Grid
            'h3grid = vgrid.grid.generator.h3grid:main',
            's2grid = vgrid.grid.generator.s2grid:main',
            'rhealpixgrid = vgrid.grid.generator.rhealpixgrid:main',
            'eaggrisea4tgrid = vgrid.grid.generator.eaggrisea4tgrid:main',
            # 'olcgrid = vgrid.grid.generator.olcgrid:main',
            'geohashgrid = vgrid.grid.generator.geohashgrid:main',           
            'gzd = vgrid.grid.generator.gzd:main',  
            'mgrsgrid = vgrid.grid.generator.mgrsgrid:main',
            'tilegrid = vgrid.grid.generator.tilegrid:main', 
            'maidenheadgrid = vgrid.grid.generator.maidenheadgrid:main',           
           
            # Grid Stats
            'h3stats = vgrid.grid.stats.h3stats:main',
            's2stats = vgrid.grid.stats.s2stats:main',
            'rhealpixstats = vgrid.grid.stats.rhealpixstats:main',
            'eaggrisea4tstats = vgrid.grid.stats.eaggrisea4tstats:main',

            # 'olcstats = vgrid.stats.olcstats:main',
            'geohashstats = vgrid.grid.stats.geohashstats:main',
            'georefstats = vgrid.grid.stats.georefstats:main',
            'mgrsstats = vgrid.grid.stats.mgrsstats:main',
            # 'tilestats = vgrid.stats.tilestats:main',
            
            'maidenheadstats = vgrid.grid.stats.maidenheadstats:main',
            'garsstats = vgrid.grid.stats.garsstats:main',

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

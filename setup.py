# python setup.py sdist bdist_wheel
# twine upload dist/*

from setuptools import setup, find_packages

requirements = [
    'tqdm~=4.66.2',
    'shapely~=2.0.1',
    'protobuf~=5.26.1',
    'fiona~=1.10.0'
],

setup(
    name='vgrid',
    version='1.0.1',
    author = 'Thang Quach',
    author_email= 'quachdongthang@gmail.com',
    url='https://github.com/thangqd/vgrid',
    description='Vgrid - A Global Geocoding System based on Vector Tiles',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    requires_python=">=3.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [            
            'vcode2geojson = vgrid.vcode:vcode2geojson_cli',  
            'vgrid = vgrid.vcode.vgrid:main',      
    
            'olcgrid = vgrid.utils.grid.olcgrid:main',
            'geohashgrid = vgrid.utils.grid.geohashgrid:main',
            'h3grid = vgrid.utils.grid.h3grid:main',
            's2grid = vgrid.utils.grid.s2grid:main',
            'maidenheadgrid = vgrid.utils.grid.maidenheadgrid:main',
            'mgrsgrid = vgrid.utils.grid.mgrsgrid:main',     
            'gzd = vgrid.utils.grid.gzd:main'            
       
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

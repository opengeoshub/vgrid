# python setup.py sdist bdist_wheel
# twine upload dist/*

from setuptools import setup, find_packages

requirements = [
    'tqdm~=4.66.2',
    'shapely~=2.0.1',
    'protobuf~=5.26.1'
],

setup(
    name='vcode',
    version='1.0.0',
    author = 'Thang Quach',
    author_email= 'quachdongthang@gmail.com',
    url='https://github.com/thangqd/vcode',
    description='Vcode - A Global Geocoding System based on Vector Tiles',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    requires_python=">=3.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [            
            'vcode2geojson = tiles_util.vcode:vcode2geojson_cli',  
            'vgrid = tiles_util.vcode.vgrid:main',      
    
            'vectortilegrid = tiles_util.utils.grid.vectortilegrid:main',
            'pluscodegrid = tiles_util.utils.grid.pluscodegrid:main',
            'geohashgrid = tiles_util.utils.grid.geohashgrid:main',
            'h3grid = tiles_util.utils.grid.h3grid:main',
            's2grid = tiles_util.utils.grid.s2grid:main',
            'maidenheadgrid = tiles_util.utils.grid.maidenheadgrid:main',
            'mgrsgrid = tiles_util.utils.grid.mgrsgrid:main'            
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

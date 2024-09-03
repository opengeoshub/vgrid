# Vgrid - A Global Geocoding System based on Vector Tiles

## Installation: 
- Using pip install (Windows/ Linux):
    ``` bash 
    pip install vgrid
    ```
- Show information of installed mbtiles-util: 
    ``` bash 
    pip show vgrid
    ```
- Install the latest vertion of mbtiles-util:
    ``` bash 
    pip install vgrid --upgrade
    ```
    
- Visit vgrid on [PyPI](https://pypi.org/project/vgrid/)

## Usage:
### vgrid:
- Create a debug Vgrid:  
    ``` bash 
    > vgrid -minzoom <minzoom>  -maxzoom <maxzoom>  -o <output_file> 
    ```
  Ex: `> vgrid -minzoom 0  -maxzoom 6>  -o vgrid0_6.mbtiles`
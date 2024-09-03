# Vcode - A Global Geocoding System based on Vector Tiles

## Installation: 
- Using pip install (Windows/ Linux):
    ``` bash 
    pip install vcode
    ```
- Show information of installed mbtiles-util: 
    ``` bash 
    pip show vcode
    ```
- Install the latest vertion of mbtiles-util:
    ``` bash 
    pip install vcode --upgrade
    ```
    
- Visit vcode on [PyPI](https://pypi.org/project/vcode/)

## Usage:
### vgrid:
- Create a debug Vgrid for Vcode:  
    ``` bash 
    > vgrid -minzoom <minzoom>  -maxzoom <maxzoom>  -o <output_file> 
    ```
  Ex: `> vgrid -minzoom 0  -maxzoom 6>  -o vgrid0_6.mbtiles`
# Vector Tiles Utilities
# Repository Renamed
**Note:** This repository has been renamed from `mbtiles-util` to `tiles-util`. Please update your references accordingly.

## Installation: 
- Using pip install (Windows/ Linux):
    ``` bash 
    pip install tiles-util
    ```
- Show information of installed mbtiles-util: 
    ``` bash 
    pip show tiles-util
    ```
- Install the latest vertion of mbtiles-util:
    ``` bash 
    pip install tiles-util --upgrade
    ```
    
- Visit tiles-util on [PyPI](https://pypi.org/project/tiles-util/)

## Usage:
### mbtilesinfo:
- Display MBTiles metadata info:  
    ``` bash 
    > mbtilesinfo <file_path>
    ```
  Ex: `> mbtilesinfo tiles.mbtiles`

### mbtiles2folder: 
- Convert MBTiles file to folder: (support raster MBTiles (.png, .jpg, .webp) and vector MBTiles (.pbf)) 
    ``` bash 
    > mbtiles2folder -i <file_name.mbtiles> -o [output_folder (optional)] -tms [TMS scheme (optional 0 or 1, default is 0)]
    ```
  Ex: `> mbtiles2folder -i tiles.mbtiles -o tiles_folder -tms 0`
  
  Without -o parameter: output_folder has the same name with MBTiles file name at current directory 

### folder2mbtiles: 
- Convert a tiles folder to MBTiles file: (support raster tile (.png, .jpg, .webp) and vector tile (.pbf))
    ``` bash 
    > folder2mbtiles -i <input_folder> -o [file_name.mbtiles (optional)] -tms [TMS scheme (optional 0 or 1, default is 0)]
    ```
  Ex: `> folder2mbtiles -i tiles_folder -o tiles.mbtiles -tms 0`
  
  Without -o parameter: file_name.mbtiles has the same name with input_folder name at current directory 

### geoson2mbtiles: 
- Convert geojson file to mbtiles (python wrapper for tippecanoe in Linux)
    ``` bash 
    > geojson2mbtiles -i <input files> -z <maxzoom> -o <output> -t <tippecanoe path> --extra-args <drop-densest-as-needed or coalesce-densest-as-needed or extend-zooms-if-still-dropping>
    ```
  Ex: ` > geojson2mbtiles -i state.geojson district.geojson -z9 -o state_district.mbtiles -t /usr/local/bin/ --extra-args coalesce-densest-as-needed`  

### folder2s3: 
- Uplpad a vector/ raster tiles folder to Amazon S3 Bucket:  
    ``` bash 
    > folder2s3 -i <input_folder>  -format <'pbf' or 'png' or 'jpg' or 'jpeg' or 'webp'>
    ```
  Ex: `> folder2s3 -i vectortiles_folder -format pbf`
 
  Input S3 parameters:

  ```bash
    > S3 Bucket name: <Your S3 Bucket name (required)>
    > S3 Prefix: [Your S3 prefix (Optional. Press Enter to upload to the bucket root folder)]
    > AWS Access Key ID: <Your_Access_Key_ID (required)>
    > AWS Secret Access Key: <Your_Secrect_Access_Key (required)>
    > AWS Region: <AWS region (Optional. Press Enter to choose default region)>
  ```
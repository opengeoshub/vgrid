# Unified Vector2H3 Function

The `vector2h3` function provides a unified interface for converting various geospatial vector data formats to H3 grid cells. It focuses on geospatial data formats and Shapely geometry objects.

## Features

- **Multiple Input Formats**: Supports GeoJSON, Parquet, Shapefile, GeoPackage, and more
- **Remote Data**: Can read from URLs
- **Flexible Input**: Accepts file paths, DataFrames, GeoJSON dictionaries, and Shapely geometries
- **Auto-detection**: Automatically detects file types and formats
- **Custom Parameters**: Supports pandas read parameters for fine-tuning
- **CLI Interface**: Command-line tools for batch processing

## Supported Input Formats

### File Formats
- **GeoJSON** (.geojson, .json) - Vector data format
- **Parquet** (.parquet) - Columnar data format with geometry support
- **Shapefile** (.shp) - ESRI shapefile format (requires geopandas)
- **GeoPackage** (.gpkg) - OGC GeoPackage format (requires geopandas)

### Data Types
- **pandas DataFrame** - In-memory tabular data with geometry column
- **GeoDataFrame** - Geospatial tabular data (requires geopandas)
- **GeoJSON dictionary** - Vector data in GeoJSON format
- **Shapely geometry objects** - Individual or list of geometries
- **URLs** - Remote files accessible via HTTP/HTTPS

## Installation

The function requires the following dependencies:
```bash
pip install pandas geopandas shapely h3 requests
```

## Basic Usage

### Python API

```python
from vgrid.conversion.vector2dggs.vector2h3 import vector2h3

# Convert GeoJSON file
result = vector2h3('data.geojson', resolution=8)

# Convert DataFrame with geometry column
import pandas as pd
from shapely.geometry import Point
df = pd.DataFrame({
    'geometry': [Point(-74.0060, 40.7128), Point(-118.2437, 34.0522)],
    'name': ['New York', 'Los Angeles']
})
result = vector2h3(df, resolution=8)

# Convert Shapely geometries
from shapely.geometry import Point, Polygon
point = Point(-74.0060, 40.7128)
result = vector2h3(point, resolution=8)

# Convert from URL
result = vector2h3('https://example.com/data.geojson', resolution=8)
```

### Command Line Interface

```bash
# Convert GeoJSON file
python -m vgrid.conversion.vector2dggs.vector2h3 -r 8 -input data.geojson

# Convert from URL
python -m vgrid.conversion.vector2dggs.vector2h3 -r 8 -input https://example.com/data.geojson

# Use stdin
cat data.geojson | python -m vgrid.conversion.vector2dggs.vector2h3 -r 8 -input -
```

## Function Parameters

### Main Parameters
- `input_data`: File path, URL, DataFrame, GeoJSON dict, or Shapely geometry
- `resolution` (int): H3 resolution [0..15]
- `compact` (bool): Enable H3 compact mode for polygons (default: False)
- `topology` (bool): Enable H3 topology preserving mode (default: False)

### Additional Parameters
- `**kwargs`: Additional arguments passed to pandas read functions

## Examples

### DataFrame with Geometry Column
```python
import pandas as pd
from shapely.geometry import Point, Polygon

data = {
    'geometry': [
        Point(-74.0060, 40.7128),
        Polygon([(-74.1, 40.6), (-73.9, 40.6), (-73.9, 40.8), (-74.1, 40.8), (-74.1, 40.6)])
    ],
    'name': ['New York', 'Manhattan Area']
}
df = pd.DataFrame(data)
result = vector2h3(df, resolution=8)
```

### GeoDataFrame
```python
import geopandas as gpd
gdf = gpd.read_file('data.shp')
result = vector2h3(gdf, resolution=8)
```

### Shapely Geometries with Properties
```python
from shapely.geometry import Point, Polygon
geometries = [Point(-74.0060, 40.7128), Polygon([...])]
properties = [{'name': 'NYC'}, {'name': 'Area'}]
result = shapely2h3(geometries, resolution=8, properties_list=properties)
```

### Compact Mode for Polygons
```python
result = vector2h3('polygons.geojson', resolution=8, compact=True)
```

### Topology Preserving Mode
```python
result = vector2h3('data.geojson', resolution=8, topology=True)
```

## Output Format

The function returns a GeoJSON FeatureCollection with H3 grid cells:

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Polygon",
        "coordinates": [[[lon1, lat1], [lon2, lat2], ...]]
      },
      "properties": {
        "h3": "8828308281fffff",
        "resolution": 8,
        "num_edges": 6,
        "original_properties": {...}
      }
    }
  ]
}
```

## Error Handling

The function provides comprehensive error handling:

- **File not found**: Clear error message for missing files
- **Invalid format**: Automatic format detection with fallback
- **Missing geometry**: Validation of geometry columns in DataFrames
- **Network errors**: Proper handling of URL access issues
- **Invalid resolution**: Validation of H3 resolution range

## Performance Considerations

- **Large files**: For very large files, consider using chunking or streaming
- **Memory usage**: GeoDataFrames and large DataFrames can be memory-intensive
- **Network**: Remote files are downloaded to memory before processing
- **Resolution**: Higher H3 resolutions generate more cells and take longer to process

## CLI Options

### Basic Options
- `-r, --resolution`: H3 resolution [0..15] (required)
- `-input, --input`: Input file path, URL, or '-' for stdin (required)

### Processing Options
- `-compact`: Enable H3 compact mode for polygons
- `-topology`: Enable H3 topology preserving mode

### Output Options
- `-output, --output`: Output file path (default: auto-generated)

## Advanced Usage

### Custom File Reading
```python
# Read with specific parameters
result = vector2h3('data.geojson', resolution=8, encoding='utf-8')
```

### Batch Processing
```python
import glob

# Process multiple GeoJSON files
for geojson_file in glob.glob('data/*.geojson'):
    result = vector2h3(geojson_file, resolution=8)
    # Save result...
```

### Integration with Other Tools
```python
# Convert and save to different formats
result = vector2h3('data.geojson', resolution=8)

# Save as GeoJSON
with open('output.geojson', 'w') as f:
    json.dump(result, f)

# Convert to GeoDataFrame
import geopandas as gpd
gdf = gpd.GeoDataFrame.from_features(result['features'])
gdf.to_file('output.gpkg', driver='GPKG')
```

## Troubleshooting

### Common Issues

1. **Missing geometry column**: Ensure DataFrame has a geometry column with Shapely objects
2. **File encoding**: Use `encoding` parameter for non-UTF-8 files
3. **Large files**: Consider using chunking for very large datasets
4. **Network issues**: Check URL accessibility and network connectivity
5. **Memory errors**: Reduce file size or use streaming for large datasets

### Debug Mode
Enable verbose output to see detailed processing information:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Contributing

The unified vector2h3 function is designed to be extensible. To add support for new file formats:

1. Add format detection in `detect_file_type()`
2. Add reading logic in `read_local_file()` or `read_remote_file()`
3. Update documentation and examples

## License

This functionality is part of the vgrid package and follows the same license terms. 
#!/usr/bin/env python3
"""
GeoJSON to CSV Converter

This script converts GeoJSON files to CSV format, extracting coordinates and properties
from GeoJSON features into a tabular format.
"""

import json
import csv
import sys
import argparse
import requests
import os
from pathlib import Path
from typing import Dict, List, Any, Union
from io import StringIO
from urllib.parse import urlparse
from tqdm import tqdm
from shapely.geometry import Polygon, box, Point, LineString
from vgrid.utils.download import is_url, read_geojson_file


def extract_coordinates(geometry: Dict[str, Any]) -> Dict[str, float]:
    """
    Extract coordinates from a GeoJSON geometry object.
    
    Args:
        geometry (Dict[str, Any]): GeoJSON geometry object
        
    Returns:
        Dict[str, float]: Dictionary containing x and y coordinates
    """
    if geometry['type'] == 'Point':
        return {'x': geometry['coordinates'][0], 'y': geometry['coordinates'][1]}
    elif geometry['type'] in ['LineString', 'MultiPoint']:
        # For LineString and MultiPoint, use the first point
        return {'x': geometry['coordinates'][0][0], 'y': geometry['coordinates'][0][1]}
    elif geometry['type'] in ['Polygon', 'MultiLineString']:
        # For Polygon and MultiLineString, use the first point of the first ring/line
        return {'x': geometry['coordinates'][0][0][0], 'y': geometry['coordinates'][0][0][1]}
    elif geometry['type'] == 'MultiPolygon':
        # For MultiPolygon, use the first point of the first ring of the first polygon
        return {'x': geometry['coordinates'][0][0][0][0], 'y': geometry['coordinates'][0][0][0][1]}
    else:
        raise ValueError(f"Unsupported geometry type: {geometry['type']}")


def geojson2csv(geojson_data):
    """
    Convert GeoJSON data to CSV format.
    
    Args:
        geojson_data (Dict): GeoJSON data as a dictionary with 'type' and 'features' keys
        
    Returns:
        str: CSV data as a string
        
    Raises:
        ValueError: If the GeoJSON format is invalid or contains no features
    """   
    features = geojson_data['features']
    if not features:
        raise ValueError("No features found in GeoJSON data")
    
    # Prepare CSV headers
    headers = ['x', 'y']  # Always include coordinates
    if features[0].get('properties'):
        headers.extend(features[0]['properties'].keys())
    
    # Create CSV in memory
    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=headers)
    writer.writeheader()
    
    for feature in features:
        row = extract_coordinates(feature['geometry'])
        if feature.get('properties'):
            row.update(feature['properties'])
        writer.writerow(row)
    
    return output.getvalue()


def geojson2csv_cli() -> None:
    """
    Command-line interface for converting GeoJSON to CSV.
    
    Usage:
        python geojson2csv.py -geojson input.geojson
        python geojson2csv.py -geojson https://example.com/data.geojson
    """
    parser = argparse.ArgumentParser(description='Convert GeoJSON to CSV format')
    parser.add_argument('-geojson', required=True, help='Input GeoJSON file path or URL')
    
    args = parser.parse_args()
    
    # Read GeoJSON data from file or URL
    geojson_data = read_geojson_file(args.geojson)
    if geojson_data is None:
        sys.exit(1)
    
    try:
        csv_data = geojson2csv(geojson_data)
        
        # Create output CSV path
        if is_url(args.geojson):
            # For URLs, use the last part of the URL as the base filename
            base_name = os.path.splitext(os.path.basename(urlparse(args.geojson).path))[0]
            if not base_name:  # If URL doesn't end with a filename
                base_name = "geojson_data"
            output_path = f"{base_name}.csv"
        else:
            # For local files, replace .geojson extension with .csv
            input_path = Path(args.geojson)
            output_path = input_path.with_suffix('.csv')
        
        # Write CSV data to file
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            f.write(csv_data)
        print(f"Successfully converted to {output_path}")
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    geojson2csv_cli() 
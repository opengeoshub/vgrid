import h3
import argparse
import os
import json
import pandas as pd
from tqdm import tqdm
from shapely.geometry import Polygon, mapping
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from pyproj import Geod

geod = Geod(ellps="WGS84")
chunk_size = 10_000  # Process by chunks of 10,000 rows

def h3_to_geojson(h3_id):
    """Convert H3 cell ID to a GeoJSON Polygon."""
    cell_boundary = h3.cell_to_boundary(h3_id)
    if cell_boundary:
        filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
        # Reverse lat/lon to lon/lat for GeoJSON compatibility
        reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
        cell_polygon = Polygon(reversed_boundary)
        
        center_lat, center_lon = h3.cell_to_latlng(h3_id)
        center_lat = round(center_lat, 7)
        center_lon = round(center_lon, 7)

        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]), 3)  # Area in square meters     
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters  
        avg_edge_len = round(cell_perimeter / 6, 3)
        if h3.is_pentagon(h3_id):
            avg_edge_len = round(cell_perimeter / 5, 3)   
        resolution = h3.get_resolution(h3_id)        
        
        return {
            "geometry": mapping(cell_polygon),
            "properties": {
                "h3": h3_id,
                "resolution": resolution,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "avg_edge_len": avg_edge_len,
                "cell_area": cell_area
            }
        }

def csv_to_h3(csv_file):
    """Convert CSV with H3 column to GeoJSON, preserving all rows."""
    output_file = os.path.join(os.getcwd(), f"{os.path.splitext(os.path.basename(csv_file))[0]}_csv2h3.geojson")
    
    features = []
    
    for chunk in pd.read_csv(csv_file, dtype={"h3": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc="Processing rows"):
            try:
                h3_id = row["h3"]
                h3_feature = h3_to_geojson(h3_id)
                if h3_feature:
                    h3_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    features.append({"type": "Feature", **h3_feature})
            except Exception as e:
                print(f"Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": features}
    with open(output_file, "w") as f:
        json.dump(geojson, f, indent=2)
    
    print(f"GeoJSON saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert CSV with H3 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'h3' column")
    args = parser.parse_args()
    
    if not os.path.exists(args.csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    csv_to_h3(args.csv)

if __name__ == "__main__":
    main()

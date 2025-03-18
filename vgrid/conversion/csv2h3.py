import os, argparse,json
import pandas as pd
from tqdm import tqdm
import h3

from shapely.geometry import Polygon, mapping
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import  chunk_size, geodesic_dggs_to_feature

def h32feature(h3_id):
    """Convert H3 cell ID to a GeoJSON Polygon."""
    cell_boundary = h3.cell_to_boundary(h3_id)
    if cell_boundary:
        filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
        # Reverse lat/lon to lon/lat for GeoJSON compatibility
        reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
        cell_polygon = Polygon(reversed_boundary)
        resolution = h3.get_resolution(h3_id)             
        num_edges = 6
        if h3.is_pentagon(h3_id):
            num_edges = 6         
        h3_feature = geodesic_dggs_to_feature("h3",h3_id,resolution,cell_polygon,num_edges)              
        
        return h3_feature
        
def main():
    parser = argparse.ArgumentParser(description="Convert CSV with H3 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'h3' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"h3": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                h3_id = row["h3"]
                h3_feature = h32feature(h3_id)
                if h3_feature:
                    h3_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(h3_feature)
            except Exception as e:
                print(f"Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2h3.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
    print(f"GeoJSON saved to {geojson_path}")

if __name__ == "__main__":
    main()

from vgrid.utils import s2
import argparse
import os
import json
import pandas as pd
from tqdm import tqdm
from shapely.geometry import Polygon, mapping
from pyproj import Geod
from vgrid.utils.antimeridian import fix_polygon

geod = Geod(ellps="WGS84")
chunk_size = 10_000  # Process by chunks of 10,000 rows


def s2_to_geojson(s2_token):
    # Create an S2 cell from the given cell ID
    cell_id = s2.CellId.from_token(s2_token)
    cell = s2.Cell(cell_id)
    if cell:
        # Get the vertices of the cell (4 vertices for a rectangular cell)
        vertices = [cell.get_vertex(i) for i in range(4)]
        # Prepare vertices in (longitude, latitude) format for Shapely
        shapely_vertices = []
        for vertex in vertices:
            lat_lng = s2.LatLng.from_point(vertex)  # Convert Point to LatLng
            longitude = lat_lng.lng().degrees  # Access longitude in degrees
            latitude = lat_lng.lat().degrees   # Access latitude in degrees
            shapely_vertices.append((longitude, latitude))

        # Close the polygon by adding the first vertex again
        shapely_vertices.append(shapely_vertices[0])  # Closing the polygon
        # Create a Shapely Polygon
        cell_polygon = fix_polygon(Polygon(shapely_vertices)) # Fix antimeridian
        lat_lng = cell_id.to_lat_lng()            
        # Extract latitude and longitude in degrees
        center_lat = round(lat_lng.lat().degrees,7)
        center_lon = round(lat_lng.lng().degrees,7)

        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)  # Area in square meters     
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters  
        avg_edge_len = round(cell_perimeter/4,3)

        return {
            "geometry": mapping(cell_polygon),
            "properties": {
                "s2": s2_token,
                "resolution": cell_id.level(),
                "center_lat": center_lat,
                "center_lon": center_lon,
                "avg_edge_len": avg_edge_len,
                "cell_area": cell_area
            }
        }


def csv_to_s2(csv_file):
    """Convert CSV with S2 column to GeoJSON"""
    output_file = os.path.join(os.getcwd(), f"{os.path.splitext(os.path.basename(csv_file))[0]}_csv2s2.geojson")
    
    features = []
    
    for chunk in pd.read_csv(csv_file, dtype={"s2": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc="Processing rows"):
            try:
                s2_id = row["s2"]
                s2_feature = s2_to_geojson(s2_id)
                if s2_feature:
                    s2_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    features.append({"type": "Feature", **s2_feature})
            except Exception as e:
                print(f"Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": features}
    with open(output_file, "w") as f:
        json.dump(geojson, f, indent=2)
    
    print(f"GeoJSON saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert CSV with H3 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 's2' column")
    args = parser.parse_args()
    
    if not os.path.exists(args.csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    csv_to_s2(args.csv)

if __name__ == "__main__":
    main()

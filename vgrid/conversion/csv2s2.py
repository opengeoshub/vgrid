import os,argparse,json
import pandas as pd
from tqdm import tqdm
from vgrid.generator.settings import chunk_size
from vgrid.utils import s2
from shapely.geometry import Polygon
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.antimeridian import fix_polygon

def s22feature(s2_token):
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
        resolution = cell_id.level()            
        num_edges = 4        
        
        s2_feature = geodesic_dggs_to_feature("s2",s2_token,resolution,cell_polygon,num_edges)              
        return s2_feature
    
def main():
    parser = argparse.ArgumentParser(description="Convert CSV with S2 column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 's2' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"s2": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                s2_id = row["s2"]
                s2_feature = s22feature(s2_id)
                if s2_feature:
                    s2_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(s2_feature)
            except Exception as e:
                print(f"Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2s2.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
    print(f"GeoJSON saved to {geojson_path}")


if __name__ == "__main__":
    main()

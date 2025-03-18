#################################################################################
#  Open-Eaggr isea3h
#################################################################################
def isea3h2feature(isea3h_id):
    if (platform.system() == 'Windows'): 
        isea3h_dggs = Eaggr(Model.isea3h)
        cell_to_shape = isea3h_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(isea3h_id),ShapeStringFormat.WKT)
        cell_to_shape_fixed = loads(fix_isea3h_wkt(cell_to_shape))
        if isea3h_id.startswith('00') or isea3h_id.startswith('09') or isea3h_id.startswith('14')\
            or isea3h_id.startswith('04') or isea3h_id.startswith('19'):
            cell_to_shape_fixed = fix_isea3h_antimeridian_cells(cell_to_shape_fixed)
        
        if cell_to_shape_fixed:
            resolution = len(isea3h_id)-2
            num_edges = 3
            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
            isea3h_feature = geodesic_dggs_to_feature("isea3h",isea3h_id,resolution,cell_polygon,num_edges)   
            return isea3h_feature

def csv2isea3h():
    parser = argparse.ArgumentParser(description="Convert CSV with isea3h column to GeoJSON")
    parser.add_argument("csv", help="Input CSV file with 'isea3h' column")
    args = parser.parse_args()
    csv = args.csv
    
    if not os.path.exists(csv):
        print(f"Error: Input file {args.csv} does not exist.")
        return
    
    geojson_features = []    
    for chunk in pd.read_csv(csv, dtype={"isea3h": str}, chunksize=chunk_size):
        for _, row in tqdm(chunk.iterrows(), total=len(chunk), desc=f"Processing {len(chunk)} rows"):
            try:
                isea3h_id = row["isea3h"]
                isea3h_feature = isea3h2feature(isea3h_id)
                if isea3h_feature:
                    isea3h_feature["properties"].update(row.to_dict())  # Append all CSV data to properties
                    geojson_features.append(isea3h_feature)
            except Exception as e:
                print(f"Skipping row {row.to_dict()}: {e}")
    
    geojson = {"type": "FeatureCollection", "features": geojson_features}
    geojson_name = os.path.splitext(os.path.basename(csv))[0]
    geojson_path = f"{geojson_name}2isea3h.geojson"

    with open(geojson_path, "w") as f:
        json.dump(geojson, f, indent=2)
    
    print(f"GeoJSON saved to {geojson_path}")

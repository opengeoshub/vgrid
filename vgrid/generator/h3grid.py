# Reference: https://observablehq.com/@claude-ducharme/h3-map
# https://h3-snow.streamlit.app/

import h3
from shapely.geometry import Polygon, box
import argparse
import json
import csv
from tqdm import tqdm
from pyproj import Geod
geod = Geod(ellps="WGS84")
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union

from vgrid.generator.settings import geodesic_dggs_to_feature
max_cells = 100_000_000

def fix_h3_antimeridian_cells(hex_boundary, threshold=-128):
    if any(lon < threshold for _, lon in hex_boundary):
        # Adjust all longitudes accordingly
        return [(lat, lon - 360 if lon > 0 else lon) for lat, lon in hex_boundary]
    return hex_boundary

def generate_grid(resolution, format='geojson'):
    base_cells = h3.get_res0_cells()
    num_base_cells = len(base_cells)
    h3_features = []
    # Progress bar for base cells
    with tqdm(total=num_base_cells, desc="Processing base cells", unit=" cells") as pbar:
        for cell in base_cells:
            child_cells = h3.cell_to_children(cell, resolution)
            # Progress bar for child cells
            for child_cell in child_cells:
                # Get the boundary of the cell
                hex_boundary = h3.cell_to_boundary(child_cell)
                # Wrap and filter the boundary
                filtered_boundary = fix_h3_antimeridian_cells(hex_boundary)
                # Reverse lat/lon to lon/lat for GeoJSON compatibility
                reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
                cell_polygon = Polygon(reversed_boundary)                
                if cell_polygon.is_valid:
                    h3_id = str(child_cell)
                    num_edges = 6
                    if (h3.is_pentagon(h3_id)):
                        num_edges = 5
                    h3_feature = geodesic_dggs_to_feature("h3",h3_id,resolution,cell_polygon,num_edges)   
                    h3_features.append(h3_feature)
                    pbar.update(1)

    if format.lower() == 'csv':
        csv_rows = []
        for feature in h3_features:
            props = feature['properties']
            row = {
                'h3': props['h3']
            }
            csv_rows.append(row)
        return csv_rows
    else:
        return {
            "type": "FeatureCollection",
            "features": h3_features,
        }
    
def geodesic_buffer(polygon, distance):
    buffered_coords = []
    for lon, lat in polygon.exterior.coords:
        # Generate points around the current vertex to approximate a circle
        circle_coords = [
            geod.fwd(lon, lat, azimuth, distance)[:2]  # Forward calculation: returns (lon, lat, back_azimuth)
            for azimuth in range(0, 360, 10)  # Generate points every 10 degrees
        ]
        buffered_coords.append(circle_coords)
    
    # Flatten the list of buffered points and form a Polygon
    all_coords = [coord for circle in buffered_coords for coord in circle]
    return Polygon(all_coords).convex_hull

def generate_grid_within_bbox(resolution, bbox, format='geojson'):
    bbox_polygon = box(*bbox)  # Create a bounding box polygon
    distance = h3.average_hexagon_edge_length(resolution,unit='m')*2
    bbox_buffer = geodesic_buffer(bbox_polygon, distance)
    bbox_buffer_cells  = h3.geo_to_cells(bbox_buffer,resolution)
    total_cells = len(bbox_buffer_cells)
    print(f"Resolution {resolution} within bounding box {bbox} will generate {total_cells} cells ")
    
    if total_cells > max_cells:
        print(f"which exceeds the limit of {max_cells}. ")
        print("Please select a smaller resolution and try again.")
        return
    else:    
        h3_features = []
        # Progress bar for base cells
        for bbox_buffer_cell in tqdm(bbox_buffer_cells, desc="Processing cells"):
                # Get the boundary of the cell
                hex_boundary = h3.cell_to_boundary(bbox_buffer_cell)
                # Wrap and filter the boundary
                filtered_boundary = fix_h3_antimeridian_cells(hex_boundary)
                # Reverse lat/lon to lon/lat for GeoJSON compatibility
                reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
                cell_polygon = Polygon(reversed_boundary)
                if cell_polygon.intersects(bbox_polygon):
                    h3_id = str(bbox_buffer_cell)
                    num_edges = 6
                    if (h3.is_pentagon(h3_id)):
                        num_edges = 5
                    h3_feature = geodesic_dggs_to_feature("h3",h3_id,resolution,cell_polygon,num_edges)   
                    h3_features.append(h3_feature)

        if format.lower() == 'csv':
            csv_rows = []
            for feature in h3_features:
                props = feature['properties']
                row = {
                    'h3': props['h3']                    
                }
                csv_rows.append(row)
            return csv_rows
        else:
            return {
                "type": "FeatureCollection",
                "features": h3_features,
            }

def generate_grid_resample(resolution, geojson_features):
    # Create a unified geometry from all features
    geometries = [shape(feature["geometry"]) for feature in geojson_features["features"]]

    unified_geom = unary_union(geometries)

    # Estimate buffer distance based on resolution
    distance = h3.average_hexagon_edge_length(resolution, unit='m') * 2
    buffered_geom = geodesic_buffer(unified_geom, distance)

    # Generate H3 cells that cover the buffered geometry
    h3_cells = h3.geo_to_cells(buffered_geom, resolution)

    h3_features = []
    for h3_cell in tqdm(h3_cells, desc="Generating H3 DGGS", unit= " cells"):
        hex_boundary = h3.cell_to_boundary(h3_cell)
        filtered_boundary = fix_h3_antimeridian_cells(hex_boundary)
        reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
        cell_polygon = Polygon(reversed_boundary)

        # Only keep cells that intersect the unified input geometry
        if cell_polygon.intersects(unified_geom):
            h3_id = str(h3_cell)
            num_edges = 6 if not h3.is_pentagon(h3_id) else 5
            h3_feature = geodesic_dggs_to_feature("h3", h3_id, resolution, cell_polygon, num_edges)
            h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }

def h3grid(resolution, bbox=None, format='geojson'):
    """
    Generate H3 grid for pure Python usage.
    
    Args:
        resolution (int): H3 resolution [0..15]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        format (str, optional): Output format ('geojson' or 'csv'). Defaults to 'geojson'.
    
    Returns:
        dict or list: GeoJSON FeatureCollection or list of CSV rows depending on format
    """
    if resolution < 0 or resolution > 15:
        raise ValueError("Resolution must be in range [0..15]")

    if bbox is None:
        bbox = [-180, -90, 180, 90]
        num_cells = h3.get_num_cells(resolution)
        if num_cells > max_cells:
            raise ValueError(f"Resolution {resolution} will generate {num_cells} cells which exceeds the limit of {max_cells}")
        return generate_grid(resolution, format)
    else:
        return generate_grid_within_bbox(resolution, bbox, format)

def h3grid_cli():
    """CLI interface for generating H3 grid."""
    parser = argparse.ArgumentParser(description="Generate H3 DGGS.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..15]")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )
    parser.add_argument(
        '-f', '--format', type=str, choices=['geojson', 'csv'], default='geojson',
        help="Output format (geojson or csv)"
    )
    args = parser.parse_args()
    
    try:
        result = h3grid(args.resolution, args.bbox, args.format)
        if result is None:
            return

        # Define the output file path
        output_path = f"h3_grid_{args.resolution}.{args.format}"
        
        if args.format == 'csv':
            with open(output_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['h3_id', 'resolution', 'num_edges', 'geometry'])
                writer.writeheader()
                writer.writerows(result)
        else:
            with open(output_path, 'w') as f:
                json.dump(result, f, indent=2)

        print(f"Output saved as {output_path}")
        
    except ValueError as e:
        print(f"Error: {str(e)}")
        return
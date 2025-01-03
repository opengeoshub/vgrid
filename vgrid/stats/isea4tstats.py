import locale
import argparse
import csv
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.enums.model import Model
from shapely.wkt import loads
from vgrid.conversion.latlon2cell import latlon2isea4t
from texttable import Texttable

from pyproj import Geod
geod = Geod(ellps="WGS84")

current_locale = locale.getlocale()  # Get the current locale setting
locale.setlocale(locale.LC_ALL, current_locale)  # Set locale to current to format numbers

isea4t_dggs = Eaggr(Model.ISEA4T)

def fix_eaggr_wkt(eaggr_wkt):
    # Extract the coordinate section
    coords_section = eaggr_wkt[eaggr_wkt.index("((") + 2 : eaggr_wkt.index("))")]
    coords = coords_section.split(",")
    # Append the first point to the end if not already closed
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    fixed_coords = ", ".join(coords)
    return f"POLYGON (({fixed_coords}))"

def isea4t_metrics(res):
    num_cells = 20*(4**res)   
    lat,lon = 10.775275567242561, 106.70679737574993
    eaggr_cell = DggsCell(latlon2isea4t(lat,lon,res))
    cell_to_shp = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell,ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_eaggr_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)
    
    isea4t_cell = DggsCell(latlon2isea4t(lat,lon,res))
    isea4t2point = isea4t_dggs.convert_dggs_cell_to_point(isea4t_cell)      
    accuracy = isea4t2point.get_accuracy()
    
    avg_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
    avg_edge_length = abs(geod.geometry_area_perimeter(cell_polygon)[1])/3  # Perimeter in meters/ 3      
    return num_cells, avg_edge_length, avg_area, accuracy


def isea4t_stats(min_res=0, max_res=39, output_file=None):
    
    t = Texttable()
    
    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)", "Accuracy"])
    
    # Check if an output file is specified (for CSV export)
    if output_file:
        # Open the output CSV file for writing
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)","Accuracy"])
            
            # Iterate through resolutions and write rows to the CSV file
            for res in range(min_res, max_res):
                num_cells, avg_edge_length, avg_area, accuracy = isea4t_metrics(res)              
                writer.writerow([res, num_cells, avg_edge_length, avg_area,accuracy])
        print(f'OpenEAGGGR ISEA4T stats saved to {output_file}.')
    else:
        for res in range(min_res, max_res+1):
            num_cells, avg_edge_length, avg_area,accuracy = isea4t_metrics(res)  
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            formatted_edge_length = locale.format_string("%.3f", avg_edge_length, grouping=True)            
            formatted_area = locale.format_string("%.3f", avg_area, grouping=True)   
            formatted_accuracy = locale.format_string("%.3f", accuracy, grouping=True)
        
            t.add_row([res, formatted_num_cells, formatted_edge_length, formatted_area,formatted_accuracy])
        
        # Print the formatted table to the console
        print(t.draw())

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display OpenEAGGR ISEA4T DGGS stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres','--minres', type=int, default=0, help="Minimum resolution.")
    parser.add_argument('-maxres','--maxres', type=int, default=39, help="Maximum resolution.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    isea4t_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()
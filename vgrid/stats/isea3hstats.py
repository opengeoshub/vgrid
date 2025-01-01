import locale
import argparse
import csv
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.enums.model import Model
from shapely.wkt import loads
from pyproj import Geod
from vgrid.conversion.latlon2cell import latlon2isea3h

from texttable import Texttable

def fix_eaggr_wkt(eaggr_wkt):
        # Extract the coordinate section
        coords_section = eaggr_wkt[eaggr_wkt.index("((") + 2 : eaggr_wkt.index("))")]
        coords = coords_section.split(",")
        # Append the first point to the end if not already closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        fixed_coords = ", ".join(coords)
        return f"POLYGON (({fixed_coords}))"

def isea3h_metrics(res):    
    num_cells = 20*(7**res)
    isea3h_dggs = Eaggr(Model.ISEA3H)
    lat,lon = 10.775275567242561, 106.70679737574993
    isea3h_cell = DggsCell(latlon2isea3h(lat,lon,res))
    cell_to_shp = isea3h_dggs.convert_dggs_cell_outline_to_shape_string(isea3h_cell,ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_eaggr_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)
    geod = Geod(ellps="WGS84")
    
    avg_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
    avg_edge_length = abs(geod.geometry_area_perimeter(cell_polygon)[1])/6  # Perimeter in meters/ 6     
    return num_cells, avg_edge_length, avg_area


def isea3h_stats(min_res=0, max_res=18, output_file=None):    
    t = Texttable()    
    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
    
    # Check if an output file is specified (for CSV export)
    if output_file:
        # Open the output CSV file for writing
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
            
            # Iterate through resolutions and write rows to the CSV file
            for res in range(min_res, max_res + 1):
                num_cells, avg_edge_length, avg_area = isea3h_metrics(res)              
                # Write to CSV without formatting locale
                writer.writerow([res, num_cells, avg_edge_length, avg_area])
        print(f'OpenEAGGGR ISEA3H stats saved to {output_file}.')
    else:
        # If no output file is provided, print the result using locale formatting in Texttable
        current_locale = locale.getlocale()  # Get the current locale setting
        locale.setlocale(locale.LC_ALL, current_locale)  # Set locale to current to format numbers
        
        # Iterate through resolutions and add rows to the table
        for res in range(min_res, max_res + 1):
            num_cells, avg_edge_length, avg_area = isea3h_metrics(res)  

            formatted_cells = locale.format_string("%d", num_cells, grouping=True)
            
            avg_edge_length = round(avg_edge_length,2)
            formatted_edge_length = locale.format_string("%.2f", avg_edge_length, grouping=True)
            
            avg_area = round(avg_area,2)
            formatted_area = locale.format_string("%.2f", avg_area, grouping=True)
            
            # Add a row to the table
            t.add_row([res, formatted_cells, formatted_edge_length, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display OpenEAGGR ISEA3H DGGS stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres','--minres', type=int, default=0, help="Minimum resolution.")
    parser.add_argument('-maxres','--maxres', type=int, default=18, help="Maximum resolution.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    isea3h_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()

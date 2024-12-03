import locale
import argparse
import csv
import math
from vgrid.geocode.s2 import CellId, LatLng, Cell
from texttable import Texttable

def avg_cell_area(res):
    """
    Calculate the average area of S2 cells at a specific level.
    
    :param level: The S2 level (0-30) for which to calculate the average area.
    :return: Average cell area in square kilometers.
    """
    # https://www.jpz.se/Html_filer/wgs_84.html
    # Total surface area of a sphere (Earth) in square kilometers
    earth_surface_area_km2 = 510_065_621.724 # 510.1 million square kilometers
        # Total number of cells at the given level
    num_cells = 6 * (4**res)

    # Average area per cell at this level
    avg_area = (earth_surface_area_km2 / num_cells)*(10**6)
    return avg_area

def avg_edge_length(res):
    """
    Calculate the average edge length of S2 cells at a specific level.
    
    :param level: The S2 level (0-30) for which to calculate the average edge length.
    :return: Average edge length in kilometers.
    """
    # Earth's radius in kilometers
    earth_radius_km = 6371.0071809 
    # https://www.jpz.se/Html_filer/wgs_84.html

    # Total number of edges at this level (4 edges per cell)
    num_edges = 6 * (4**res) * 4

    # Total arc length of Earth's circumference
    earth_circumference_km = 2 * math.pi * earth_radius_km
    # earth_circumference_km = (40075.017 + 40007.863)/2

    # Average edge length for a cell at this level
    avg_edge_length = (earth_circumference_km / math.sqrt(6* (4**res)))*1000 # in meters
    return avg_edge_length


def s2_stats(min_res=0, max_res=30, output_file=None):

    # Create a Texttable object for displaying in the terminal
    t = Texttable()
    
    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(["Resolution", "Number of Cells",  "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
    
    # Check if an output file is specified (for CSV export)
    if output_file:
        # Open the output CSV file for writing
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
            
            # Iterate through resolutions and write rows to the CSV file
            for res in range(min_res, max_res + 1):
                num_cells_at_res =  6 * (4**res)
                cell_area = round(avg_cell_area(res),2)                
                cell_width = round(math.sqrt(cell_area),2)
                # Write to CSV without formatting locale
                writer.writerow([res, num_cells_at_res, cell_width, cell_area])
    else:
        # If no output file is provided, print the result using locale formatting in Texttable
        current_locale = locale.getlocale()  # Get the current locale setting
        locale.setlocale(locale.LC_ALL, current_locale)  # Set locale to current to format numbers
        
        # Iterate through resolutions and add rows to the table
        for res in range(min_res, max_res + 1):
            num_cells_at_res = 6 * (4**res)
            formatted_cells = locale.format_string("%d", num_cells_at_res, grouping=True)
            
            cell_area = round(avg_cell_area(res),2)
            formatted_area = locale.format_string("%.2f", cell_area, grouping=True)
            
            cell_width = round(math.sqrt(cell_area),2)
            formatted_width = locale.format_string("%.2f", cell_width, grouping=True)

            # Add a row to the table
            t.add_row([res, formatted_cells, formatted_width, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display S2 DGGS stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres','--minres', type=int, default=0, help="Minimum resolution.")
    parser.add_argument('-maxres','--maxres', type=int, default=22, help="Maximum resolution.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    s2_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()

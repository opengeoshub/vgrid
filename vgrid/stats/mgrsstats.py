import locale
import argparse
import csv
import math
from texttable import Texttable

def mgrs_metrics(res):
    latitude_degrees = 8  # The latitude span of each GZD cell in degrees
    longitude_degrees = 6  # The longitude span of each GZD cell in degrees
    km_per_degree = 111 # Approximate kilometers per degree of latitude/longitude
    gzd_cells = 1200  # Total number of GZD cells

    # Convert degrees to kilometers
    latitude_span = latitude_degrees * km_per_degree
    longitude_span = longitude_degrees * km_per_degree

    # Map resolution levels to corresponding cell sizes (in kilometers)
    res_to_cell_size = {
        0: 100,  # 100 km x 100 km
        1: 10,   # 10 km x 10 km
        2: 1,    # 1 km x 1 km
        3: 0.1,  # 0.1 km x 0.1 km
        4: 0.01, # 0.01 km x 0.01 km
        5: 0.001 # 0.001 km x 0.001 km
    }

    cell_size = res_to_cell_size[res]
    # Calculate number of cells in latitude and longitude for the chosen cell size
    cells_latitude = latitude_span / cell_size
    cells_longitude = longitude_span / cell_size

    # Total number of cells for each GZD cell
    cells_per_gzd_cell = cells_latitude * cells_longitude

    # Total number of cells for all GZD cells
    num_cells = cells_per_gzd_cell * gzd_cells
    avg_edge_length = cell_size*(10**3)
    avg_area = avg_edge_length**2
       
    return num_cells, avg_edge_length, avg_area


def mgrs_stats(min_res=0, max_res=5, output_file=None):
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
                num_cells, avg_edge_length, avg_area = mgrs_metrics(res)
                avg_area = round(avg_area,2)    
                # Write to CSV without formatting locale
                writer.writerow([res, num_cells, avg_edge_length, avg_area])
    else:
        # If no output file is provided, print the result using locale formatting in Texttable
        current_locale = locale.getlocale()  # Get the current locale setting
        locale.setlocale(locale.LC_ALL, current_locale)  # Set locale to current to format numbers
        
        # Iterate through resolutions and add rows to the table
        for res in range(min_res, max_res + 1):
            num_cells, avg_edge_length, avg_area = mgrs_metrics(res)
           
            formatted_cells = locale.format_string("%d", num_cells, grouping=True)
            
            formatted_length = locale.format_string("%.2f", avg_edge_length, grouping=True)

            avg_area = round(avg_area,2)    
            formatted_area = locale.format_string("%.2f", avg_area, grouping=True)
            
            # Add a row to the table
            t.add_row([res, formatted_cells, formatted_length, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())
        

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display mgrs stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres','--minres', type=int, default=0, help="Minimum resolution.")
    parser.add_argument('-maxres','--maxres', type=int, default=5, help="Maximum resolution.")
    args = parser.parse_args()

    print('Resolution 0: 100 x 100 km')
    print('Resolution 1: 10 x 10 km')
    print('2 <= Resolution <= 5 = Finer subdivisions (1 x 1 km, 0.1 x 0.11 km, etc.)')

    # Call the function with the provided output file (if any)
    mgrs_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()

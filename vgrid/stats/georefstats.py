import locale
import argparse
import csv
import math
from texttable import Texttable

locale.setlocale(locale.LC_ALL, '')

def georef_metrics(res):
    earth_surface_area_km2 = 510_065_621.724 
    if res == 0:
        num_cells =  (180/15)*(360/15) # 15 x 15 degrees
    elif res == 1:
        num_cells =  180*360   # 1 x 1 degrees
    elif res == 2:
        num_cells = (180*60)*(360 * 60)  # Subdivision into 1' x 1' cells
    elif res ==3:
        num_cells =  (180*600)*(360 * 600) # Subdivision into 0.1' x 0.1' cells
    elif res ==4:
        num_cells =  (180*6000)*(360 * 6000) # Subdivision into 0.01' x 0.01' cells
  
    avg_area = (earth_surface_area_km2 / num_cells)*(10**6)
    avg_edge_length = math.sqrt(avg_area)
    return num_cells, avg_edge_length, avg_area


def georef_stats(output_file=None):
    min_res=0 
    max_res=4 
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
            
            for res in range(min_res, max_res + 1):
                num_cells, avg_edge_length, avg_area = georef_metrics(res)
                writer.writerow([res, num_cells, avg_edge_length, avg_area])
    else:
        for res in range(min_res, max_res + 1):
            num_cells, avg_edge_length, avg_area = georef_metrics(res)
           
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            formatted_length = locale.format_string("%.3f", avg_edge_length, grouping=True)
            formatted_area = locale.format_string("%.3f", avg_area, grouping=True)
            # Add a row to the table
            t.add_row([res, formatted_num_cells, formatted_length, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())
        

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display GEOREF stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    args = parser.parse_args()

    print('Resolution 0: 15 x 15 degrees')
    print('Resolution 1: 1 x 1 degree')
    print('Resolution 2: 1 x 1 minute')
    print('Resolution 3: 0.1 x 0.1 minute')
    print('Resolution 4: 0.01 x 0.01 minute')

    # Call the function with the provided output file (if any)
    georef_stats(args.output)

if __name__ == "__main__":
    main()

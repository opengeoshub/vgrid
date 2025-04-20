import locale
import argparse
import csv
import h3
from texttable import Texttable

locale.setlocale(locale.LC_ALL, '')

def h3_stats(output_file=None):
    min_res=0
    max_res=15 

    # Create a Texttable object for displaying in the terminal
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
                num_cells_at_res = h3.get_num_cells(res)
                cell_width = h3.average_hexagon_edge_length(res,unit='m')
                cell_area = h3.average_hexagon_area(res,unit='m^2')
                
                # Write to CSV without formatting locale
                writer.writerow([res, num_cells_at_res, cell_width, cell_area])
    else:
        for res in range(min_res, max_res + 1):
            num_cells = h3.get_num_cells(res)
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            
            cell_width = h3.average_hexagon_edge_length(res,unit='m')
            formatted_width = locale.format_string("%.3f", cell_width, grouping=True)
            
            cell_area = h3.average_hexagon_area(res,unit='m^2')
            formatted_area = locale.format_string("%.3f", cell_area, grouping=True)
            
            # Add a row to the table
            t.add_row([res, formatted_num_cells, formatted_width, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display H3 DGGS stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    h3_stats(args.output)

if __name__ == "__main__":
    main()

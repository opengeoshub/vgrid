import locale
import argparse
import csv
import h3
from texttable import Texttable

def h3_stats(min_res=0, max_res=15, output_file=None):

    # Create a Texttable object for displaying in the terminal
    t = Texttable()
    
    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(["Resolution", "Number of Cells", "Avg Hex Length (m)", "Avg Hex Area (sq m)"])
    
    # Check if an output file is specified (for CSV export)
    if output_file:
        # Open the output CSV file for writing
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Resolution", "Number of Cells", "Avg Hex Edge Length (m)", "Avg Hex Area (sq m)"])
            
            # Iterate through resolutions and write rows to the CSV file
            for res in range(min_res, max_res + 1):
                num_cells_at_res = h3.get_num_cells(res)
                cell_width = round(h3.average_hexagon_edge_length(res,unit='m'),2)
                cell_area = round(h3.average_hexagon_area(res,unit='m^2'),2)
                
                # Write to CSV without formatting locale
                writer.writerow([res, num_cells_at_res, cell_width, cell_area])
    else:
        # If no output file is provided, print the result using locale formatting in Texttable
        current_locale = locale.getlocale()  # Get the current locale setting
        locale.setlocale(locale.LC_ALL, current_locale)  # Set locale to current to format numbers
        
        # Iterate through resolutions and add rows to the table
        for res in range(min_res, max_res + 1):
            num_cells_at_res = h3.get_num_cells(res)
            formatted_cells = locale.format_string("%d", num_cells_at_res, grouping=True)
            
            cell_width = round(h3.average_hexagon_edge_length(res,unit='m'),2)
            formatted_width = locale.format_string("%.2f", cell_width, grouping=True)
            
            cell_area = round(h3.average_hexagon_area(res,unit='m^2'),2)
            formatted_area = locale.format_string("%.2f", cell_area, grouping=True)
            
            # Add a row to the table
            t.add_row([res, formatted_cells, formatted_width, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display H3 DGGS stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres','--minres', type=int, default=0, help="Minimum resolution.")
    parser.add_argument('-maxres','--maxres', type=int, default=15, help="Maximum resolution.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    h3_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()

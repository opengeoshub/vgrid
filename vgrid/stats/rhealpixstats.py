import locale
import argparse
import csv
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round, wrap_longitude, wrap_latitude
from texttable import Texttable

def rheapix_stats(min_res=0, max_res=15, output_file=None):
    rdggs = RHEALPixDGGS()
    
    # Create a Texttable object for displaying in the terminal
    t = Texttable()
    
    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(["Resolution", "Number of Cells", "Cell Width (m)", "Cell Area (sq m)"])
    
    # Check if an output file is specified (for CSV export)
    if output_file:
        # Open the output CSV file for writing
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Resolution", "Number of Cells", "Cell Width (m)", "Cell Area (sq m)"])
            
            # Iterate through resolutions and write rows to the CSV file
            for res in range(min_res, max_res + 1):
                num_cells_at_res = rdggs.num_cells(res)
                cell_width = round(rdggs.cell_width(res),2)
                cell_area = round(rdggs.cell_area(res),2)
                
                # Write to CSV without formatting locale
                writer.writerow([res, num_cells_at_res, cell_width, cell_area])
    else:
        # If no output file is provided, print the result using locale formatting in Texttable
        current_locale = locale.getlocale()  # Get the current locale setting
        locale.setlocale(locale.LC_ALL, current_locale)  # Set locale to current to format numbers
        
        # Iterate through resolutions and add rows to the table
        for res in range(min_res, max_res + 1):
            num_cells_at_res = rdggs.num_cells(res)
            formatted_cells = locale.format_string("%d", num_cells_at_res, grouping=True)
            
            cell_width = round(rdggs.cell_width(res),2)
            formatted_width = locale.format_string("%.2f", cell_width, grouping=True)
            
            cell_area = round(rdggs.cell_area(res),2)
            formatted_area = locale.format_string("%.2f", cell_area, grouping=True)
            
            # Add a row to the table
            t.add_row([res, formatted_cells, formatted_width, formatted_area])
        
        # Print the formatted table to the console
        print(t.draw())

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display RHEALPix DGGS stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres','--minres', type=int, default=0, help="Minimum resolution.")
    parser.add_argument('-maxres','--maxres', type=int, default=15, help="Maximum resolution.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    rheapix_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()

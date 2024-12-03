import math
import csv
import argparse
import locale
from texttable import Texttable

def olc_metrics(res):
    earth_surface_area_km2 = 510_065_621.724 # 510.1 million square kilometers
    base_cells, num_cells= 162, 162
    if res == 1:
        num_cells = base_cells
    elif 2 <= res <= 4:
        # Each cell is subdivided into 400 cells at ress 2 to 4
        num_cells =  base_cells * (400 ** (res - 2))
    elif res >= 5:
        # For res 5 and above, each cell is subdivided into 20 cells
        num_cells =  base_cells * (400 ** 3) * (20 ** (res - 5))
    
    avg_area = (earth_surface_area_km2 / num_cells)*(10**6)
    avg_edge_length = math.sqrt(avg_area)
    
    return num_cells, avg_edge_length, avg_area


# Function to display and/or save statistics
def olc_stats(min_res=1, max_res=15, output_file=None):
    # Create a Texttable object for displaying in the terminal
    t = Texttable()
    t.add_row(["res", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
    
    # Check if an output file is specified (for CSV export)
    if output_file:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["res", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
            for res in range(min_res, max_res + 1):
                total_cells, avg_edge_length, avg_area = olc_metrics(res)
                writer.writerow([res, total_cells, round(avg_edge_length, 2), round(avg_area, 2)])
    else:
        # Print table to console
        for res in range(min_res, max_res + 1):
            total_cells, avg_edge_length, avg_area = olc_metrics(res)
            formatted_cells = locale.format_string("%d", total_cells, grouping=True)
            formatted_edge_length = locale.format_string("%.2f", avg_edge_length, grouping=True)
            formatted_area = locale.format_string("%.2f", avg_area, grouping=True)
            t.add_row([res, formatted_cells, formatted_edge_length, formatted_area])
        print(t.draw())

# Main function to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description="Export or display OLC stats.")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    parser.add_argument('-minres', '--minres', type=int, default=1, help="Minimum res.")
    parser.add_argument('-maxres', '--maxres', type=int, default=15, help="Maximum res for OLC.")
    args = parser.parse_args()
    olc_stats(args.minres, args.maxres, args.output)

if __name__ == "__main__":
    main()

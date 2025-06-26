import locale
import argparse
import csv
import math
from texttable import Texttable

locale.setlocale(locale.LC_ALL, "")


def maidenhead_metrics(res):
    earth_surface_area_km2 = 510_065_621.724
    if res == 1:
        lon_width, lat_width = 20, 10
    elif res == 2:
        lon_width, lat_width = 2, 1
    elif res == 3:
        lon_width, lat_width = 0.083333, 0.041666  # ~5 min x 2.5 min
    elif res == 4:
        lon_width, lat_width = 0.008333, 0.004167  # ~30 sec x 15 sec

    num_lon_cells = int(360 / lon_width)
    num_lat_cells = int(180 / lat_width)
    num_cells = num_lon_cells * num_lat_cells

    avg_area = (earth_surface_area_km2 / num_cells) * (10**6)
    avg_edge_length = math.sqrt(avg_area)

    return num_cells, avg_edge_length, avg_area


def maidenhead_stats(output_file=None):
    min_res = 1
    max_res = 4
    # Create a Texttable object for displaying in the terminal
    t = Texttable()

    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(
        ["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"]
    )

    # Check if an output file is specified (for CSV export)
    if output_file:
        # Open the output CSV file for writing
        with open(output_file, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                [
                    "Resolution",
                    "Number of Cells",
                    "Avg Edge Length (m)",
                    "Avg Cell Area (sq m)",
                ]
            )

            # Iterate through resolutions and write rows to the CSV file
            for res in range(min_res, max_res + 1):
                num_cells, avg_edge_length, avg_area = maidenhead_metrics(res)
                writer.writerow([res, num_cells, avg_edge_length, avg_area])
    else:
        for res in range(min_res, max_res + 1):
            num_cells, avg_edge_length, avg_area = maidenhead_metrics(res)
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            formatted_edge_length = locale.format_string(
                "%.2f", avg_edge_length, grouping=True
            )
            formatted_area = locale.format_string("%.2f", avg_area, grouping=True)
            t.add_row([res, formatted_num_cells, formatted_edge_length, formatted_area])
        # Print the formatted table to the console
        print(t.draw())


def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display Maidenhead stats.")
    parser.add_argument("-o", "--output", help="Output CSV file name.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    maidenhead_stats(args.output)


if __name__ == "__main__":
    main()

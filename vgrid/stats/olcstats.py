import math
import csv
import argparse
import locale
from texttable import Texttable
from vgrid.generator.olcgrid import calculate_total_cells
from pyproj import Geod

geod = Geod(ellps="WGS84")

locale.setlocale(locale.LC_ALL, "")


def olc_metrics(res):
    earth_surface_area_km2 = 510_065_621.724  # 510.1 million square kilometers
    bbox = [-180, -90, 180, 90]
    num_cells = calculate_total_cells(res, bbox)
    avg_area = (earth_surface_area_km2 / num_cells) * (10**6)
    avg_edge_length = math.sqrt(avg_area)
    return num_cells, avg_edge_length, avg_area


def olc_stats(output_file=None):
    """
    Display and/or save OLC statistics to a CSV file for specific resolutions.
    """
    import locale

    locale.setlocale(locale.LC_ALL, "")  # Set locale for formatting numbers
    t = Texttable()
    t.add_row(
        ["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"]
    )

    # List of specific resolutions
    resolutions = [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]

    # Prepare CSV output if specified
    if output_file:
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
            for res in resolutions:
                num_cells, avg_edge_len, avg_area = olc_metrics(res)
                writer.writerow([res, num_cells, avg_edge_len, avg_area])
    else:
        # Display in terminal
        for res in resolutions:
            num_cells, avg_edge_len, avg_area = olc_metrics(res)
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            formatted_avg_edge_len = locale.format_string(
                "%.3f", avg_edge_len, grouping=True
            )
            formatted_avg_area = locale.format_string("%.3f", avg_area, grouping=True)
            t.add_row(
                [res, formatted_num_cells, formatted_avg_edge_len, formatted_avg_area]
            )
        print(t.draw())


def main():
    """
    Main function to handle command-line arguments and invoke OLC stats generation.
    """
    parser = argparse.ArgumentParser(
        description="Generate statistics for Open Location Code (OLC)."
    )
    parser.add_argument("-o", "--output", help="Output CSV file name.")
    args = parser.parse_args()

    # Generate stats
    olc_stats(args.output)


if __name__ == "__main__":
    main()

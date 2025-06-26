import locale
import argparse
import csv
from vgrid.utils.easedggs.constants import levels_specs
from texttable import Texttable

locale.setlocale(locale.LC_ALL, "")


def ease_stats(output_file=None):
    min_res = 0
    max_res = 6

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
                if res in levels_specs:
                    num_cells_at_res = (
                        levels_specs[res]["n_row"] * levels_specs[res]["n_col"]
                    )
                    cell_width = levels_specs[res]["x_length"]
                    cell_area = (
                        levels_specs[res]["x_length"] * levels_specs[res]["y_length"]
                    )
                    # Write to CSV without formatting locale
                    writer.writerow([res, num_cells_at_res, cell_width, cell_area])
    else:
        for res in range(min_res, max_res + 1):
            if res in levels_specs:
                num_cells = levels_specs[res]["n_row"] * levels_specs[res]["n_col"]
                formatted_num_cells = locale.format_string(
                    "%d", num_cells, grouping=True
                )

                cell_width = levels_specs[res]["x_length"]
                formatted_width = locale.format_string(
                    "%.3f", cell_width, grouping=True
                )

                cell_area = (
                    levels_specs[res]["x_length"] * levels_specs[res]["y_length"]
                )
                formatted_area = locale.format_string("%.3f", cell_area, grouping=True)

                # Add a row to the table
                t.add_row([res, formatted_num_cells, formatted_width, formatted_area])

        # Print the formatted table to the console
        print(t.draw())


def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Export or display EASE-DGGS stats.")
    parser.add_argument("-o", "--output", help="Output CSV file name.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    ease_stats(args.output)


if __name__ == "__main__":
    main()

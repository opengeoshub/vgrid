import locale
import argparse
import csv
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from texttable import Texttable

locale.setlocale(locale.LC_ALL, "")

rdggs = RHEALPixDGGS()


def rhealpix_metrics(res):
    num_cells = rdggs.num_cells(res)
    avg_edge_length = rdggs.cell_width(res)
    avg_area = rdggs.cell_area(res)
    return num_cells, avg_edge_length, avg_area


def rheapix_stats(output_file=None):
    min_res = 0
    max_res = 15
    # Create a Texttable object for displaying in the terminal
    t = Texttable()

    # Add header to the table, including the new 'Cell Width' and 'Cell Area' columns
    t.add_row(
        ["Resolution", "Number of Cells", "Avg Edge Length (m)", "Cell Area (sq m)"]
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
                    "Cell Area (sq m)",
                ]
            )
            for res in range(min_res, max_res + 1):
                num_cells, avg_edge_length, avg_area = rhealpix_metrics(res)
                writer.writerow([res, num_cells, avg_edge_length, avg_area])
        print(f"RHEALPix stats saved to {output_file}.")

    else:
        for res in range(min_res, max_res + 1):
            num_cells, avg_edge_length, avg_area = rhealpix_metrics(res)
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            formatted_edge_length = locale.format_string(
                "%.3f", avg_edge_length, grouping=True
            )
            formatted_area = locale.format_string("%.3f", avg_area, grouping=True)

            t.add_row([res, formatted_num_cells, formatted_edge_length, formatted_area])

        print(t.draw())


def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Export or display RHEALPix DGGS stats."
    )
    parser.add_argument("-o", "--output", help="Output CSV file name.")
    args = parser.parse_args()

    # Call the function with the provided output file (if any)
    rheapix_stats(args.output)


if __name__ == "__main__":
    main()

import math
import csv
import argparse
import locale
from texttable import Texttable
from vgrid.utils import olc
from vgrid.conversion.latlon2dggs import latlon2olc
from vgrid.generator.olcgrid import calculate_total_cells
from shapely.geometry import Polygon,mapping
from pyproj import Geod
geod = Geod(ellps="WGS84")

locale.setlocale(locale.LC_ALL, '')


def olc_metrics(res):
    lat,lon = 10.775275567242561, 106.70679737574993
    olc_code = latlon2olc(lat,lon, res)
    coord = olc.decode(olc_code)
    # Create the bounding box coordinates for the polygon
    min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
    max_lat, max_lon = coord.latitudeHi, coord.longitudeHi
    # Define the polygon based on the bounding box
    cell_polygon = Polygon([
        [min_lon, min_lat],  # Bottom-left corner
        [max_lon, min_lat],  # Bottom-right corner
        [max_lon, max_lat],  # Top-right corner
        [min_lon, max_lat],  # Top-left corner
        [min_lon, min_lat]   # Closing the polygon (same as the first point)
    ])
    
    avg_area = abs(geod.geometry_area_perimeter(cell_polygon)[0]) # Area in square meters     
    # Calculate width (longitude difference at a constant latitude)
    cell_width = geod.line_length([min_lon, max_lon], [min_lat, min_lat])
    
    # Calculate height (latitude difference at a constant longitude)
    cell_height = geod.line_length([min_lon, min_lon], [min_lat, max_lat])
    
    avg_edge_len = round((cell_width+cell_height)/2,3)
    
    bbox = [-180, -90, 180, 90]
    num_cells = calculate_total_cells(res, bbox)
    
    return num_cells, avg_edge_len, avg_area


def olc_stats(output_file=None):
    """
    Display and/or save OLC statistics to a CSV file for specific resolutions.
    """
    import locale
    locale.setlocale(locale.LC_ALL, '')  # Set locale for formatting numbers
    t = Texttable()
    t.add_row(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])

    # List of specific resolutions
    resolutions = [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]

    # Prepare CSV output if specified
    if output_file:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Resolution", "Number of Cells", "Avg Edge Length (m)", "Avg Cell Area (sq m)"])
            for res in resolutions:
                num_cells, avg_edge_len, avg_area = olc_metrics(res)
                writer.writerow([res, num_cells, avg_edge_len, avg_area])
    else:
        # Display in terminal
        for res in resolutions:
            num_cells, avg_edge_len, avg_area  = olc_metrics(res)
            formatted_num_cells = locale.format_string("%d", num_cells, grouping=True)
            formatted_avg_edge_len = locale.format_string("%.3f", avg_edge_len, grouping=True)
            formatted_avg_area = locale.format_string("%.3f", avg_area, grouping=True)
            t.add_row([res, formatted_num_cells, formatted_avg_edge_len, formatted_avg_area])
        print(t.draw())

def main():
    """
    Main function to handle command-line arguments and invoke OLC stats generation.
    """
    parser = argparse.ArgumentParser(description="Generate statistics for Open Location Code (OLC).")
    parser.add_argument('-o', '--output', help="Output CSV file name.")
    args = parser.parse_args()

    # Generate stats
    olc_stats(args.output)

if __name__ == "__main__":
    main()

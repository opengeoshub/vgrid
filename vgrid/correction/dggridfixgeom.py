import json
import os
import argparse
from vgrid.utils.antimeridian import fix_geojson


def fix_geom(input_file, output_file):
    with open(input_file, "r") as infile:
        geojson_data = json.load(infile)  # Load as dictionary

    fixed_geojson_data = fix_geojson(geojson_data)
    with open(output_file, "w") as outfile:
        json.dump(fixed_geojson_data, outfile, indent=2)

    print(f"GeoJSON has been fixed antimeridian and saved to {output_file}")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Fix antimeridian GeoJSON created from DGGRID"
    )
    parser.add_argument("input", help="Input GeoJSON file")

    # Parse the arguments
    args = parser.parse_args()

    # Derive output filename by appending "_fixed" to the input filename
    # output_file = f"{args.input.rsplit('.', 1)[0]}_fixed.geojson"
    output_file = os.path.join(
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(args.input))[0]}_geom_fixed.geojson",
    )

    # Call the function with the input and output file
    fix_geom(args.input, output_file)


if __name__ == "__main__":
    main()

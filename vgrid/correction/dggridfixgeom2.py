import json
import os
import argparse


def fix_antimeridian_cells(hex_boundary, threshold=-128):
    if any(lon < threshold for _, lon in hex_boundary):
        # Adjust all longitudes accordingly
        return [(lat, lon - 360 if lon > 0 else lon) for lat, lon in hex_boundary]
    return hex_boundary


def fix_geojson(geojson_data, threshold=-128):
    """
    Apply the fix_antimeridian_cells function to all geometries in the GeoJSON data.
    """
    for feature in geojson_data.get("features", []):
        geometry = feature.get("geometry", {})
        if geometry["type"] == "Polygon":
            # Fix the exterior ring
            geometry["coordinates"][0] = fix_antimeridian_cells(
                geometry["coordinates"][0], threshold
            )
            # Fix any interior rings (holes)
            for i in range(1, len(geometry["coordinates"])):
                geometry["coordinates"][i] = fix_antimeridian_cells(
                    geometry["coordinates"][i], threshold
                )
        elif geometry["type"] == "MultiPolygon":
            for polygon in geometry["coordinates"]:
                # Fix the exterior ring
                polygon[0] = fix_antimeridian_cells(polygon[0], threshold)
                # Fix any interior rings (holes)
                for i in range(1, len(polygon)):
                    polygon[i] = fix_antimeridian_cells(polygon[i], threshold)
    return geojson_data


def fix_geom(input_file, output_file):
    with open(input_file, "r") as infile:
        geojson_data = json.load(infile)  # Load as dictionary

    fixed_geojson_data = fix_geojson(geojson_data)
    with open(output_file, "w") as outfile:
        json.dump(fixed_geojson_data, outfile, indent=2)

    print(f"GeoJSON has been fixed for the antimeridian and saved to {output_file}")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Fix antimeridian GeoJSON created from DGGRID"
    )
    parser.add_argument("input", help="Input GeoJSON file")

    # Parse the arguments
    args = parser.parse_args()

    # Derive output filename by appending "_fixed" to the input filename
    output_file = os.path.join(
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(args.input))[0]}_geom_fixed_2.geojson",
    )

    # Call the function with the input and output file
    fix_geom(args.input, output_file)


if __name__ == "__main__":
    main()

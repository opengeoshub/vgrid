from shapely.geometry import Polygon, MultiPolygon, mapping
import json
import os
from vgrid.utils.antimeridian import fix_polygon


def construct_geometry(coords, is_multipolygon=False):
    if is_multipolygon:
        polygon = MultiPolygon([Polygon(poly) for poly in coords])
    else:
        polygon = Polygon(coords)

    fixed_polygon = fix_polygon(polygon)
    return fixed_polygon


def main():
    out_file_dir = os.getcwd()
    out_file_name = "cube.geojson"
    out_file = os.path.join(out_file_dir, out_file_name)

    # Define cube faces with coordinates and s2_tokens
    facets = [
        {
            "zoneID": "1",
            "coordinates": [
                [-45, -35.264389682754654],
                [45, -35.264389682754654],
                [45, 35.264389682754654],
                [-45, 35.264389682754654],
                [-45, -35.264389682754654],
            ],
        },
        {
            "zoneID": "3",
            "coordinates": [
                [45, -35.264389682754654],
                [135, -35.264389682754654],
                [135, 35.264389682754654],
                [45, 35.264389682754654],
                [45, -35.264389682754654],
            ],
        },
        {
            "zoneID": "5",
            "coordinates": [
                [45, 35.264389682754654],
                [135, 35.264389682754654],
                [-135, 35.264389682754654],
                [-45, 35.264389682754654],
                [45, 35.264389682754654],
            ],
        },
        {
            "zoneID": "7",
            "coordinates": [
                [135, 35.264389682754654],
                [135, -35.264389682754654],
                [-135, -35.264389682754654],
                [-135, 35.264389682754654],
                [135, 35.264389682754654],
            ],
        },
        {
            "zoneID": "9",
            "coordinates": [
                [-135, 35.264389682754654],
                [-135, -35.264389682754654],
                [-45, -35.264389682754654],
                [-45, 35.264389682754654],
                [-135, 35.264389682754654],
            ],
        },
        {
            "zoneID": "b",
            "coordinates": [
                [-135, -35.264389682754654],
                [135, -35.264389682754654],
                [45, -35.264389682754654],
                [-45, -35.264389682754654],
                [-135, -35.264389682754654],
            ],
        },
    ]

    # Create GeoJSON features
    geojson_features = []
    for facet in facets:
        geometry = construct_geometry(facet["coordinates"])
        geojson_features.append(
            {
                "type": "Feature",
                "geometry": mapping(geometry),
                "properties": {"zoneID": facet["zoneID"]},
            }
        )

    # Final GeoJSON output
    geojson_output = {"type": "FeatureCollection", "features": geojson_features}

    # Save as GeoJSON
    with open(out_file, "w") as f:
        json.dump(geojson_output, f, indent=2)

    print(f"Cube saved as {out_file}")


if __name__ == "__main__":
    main()

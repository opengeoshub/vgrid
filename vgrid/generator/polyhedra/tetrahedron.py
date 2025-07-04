from shapely.geometry import Polygon
import json
import os


def constructGeometry(facet):
    # Create a Polygon with the vertices (longitude, latitude)
    poly = Polygon([(v[0], v[1]) for v in facet])  # (lon, lat)
    return poly


def main():
    outFileDir = os.getcwd()  # Use current directory
    outFileName = "tetrahedron.geojson"
    outFile = os.path.join(outFileDir, outFileName)

    # Define facets with coordinates and zoneIDs
    facets = [
        {
            "zoneID": "0",
            "coordinates": [
                [-180.0, 0.0],
                [-180.0, 90.0],
                [-90.0, 90.0],
                [0.0, 90.0],
                [0.0, 0.0],
                [-90.0, 0.0],
                [-180.0, 0.0],
            ],
        },
        {
            "zoneID": "1",
            "coordinates": [
                [0.0, 0.0],
                [0.0, 90.0],
                [90.0, 90.0],
                [180.0, 90.0],
                [180.0, 0.0],
                [90.0, 0.0],
                [0.0, 0.0],
            ],
        },
        {
            "zoneID": "2",
            "coordinates": [
                [-180.0, -90.0],
                [-180.0, 0.0],
                [-90.0, 0.0],
                [0.0, 0.0],
                [0.0, -90.0],
                [-90.0, -90.0],
                [-180.0, -90.0],
            ],
        },
        {
            "zoneID": "3",
            "coordinates": [
                [0.0, -90.0],
                [0.0, 0.0],
                [90.0, 0.0],
                [180.0, 0.0],
                [180.0, -90.0],
                [90.0, -90.0],
                [0.0, -90.0],
            ],
        },
    ]

    geojson_features = []  # Store GeoJSON features

    for facet in facets:
        geojson_features.append(
            {
                "type": "Feature",
                "geometry": {"type": "Polygon", "coordinates": [facet["coordinates"]]},
                "properties": {"zoneID": facet["zoneID"]},
            }
        )

    # Save as GeoJSON with CRS
    geojson_output = {
        "type": "FeatureCollection",
        "name": "tetrahedron",
        "crs": {
            "type": "name",
            "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"},
        },
        "features": geojson_features,
    }

    with open(outFile, "w") as f:
        json.dump(geojson_output, f)

    print(f"Tetrahedron saved as {outFile}")


if __name__ == "__main__":
    main()

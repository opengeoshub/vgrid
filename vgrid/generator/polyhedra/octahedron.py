# Reference: https://github.com/paulojraposo/QTM/blob/master/qtmgenerator.py

from shapely.geometry import Polygon, LinearRing
import json
import os

def constructGeometry(facet):  
    vertexTuples = facet[:4]
    # Create a LinearRing with the vertices
    ring = LinearRing([(vT[1], vT[0]) for vT in vertexTuples])  # sequence: lon, lat (x,y)
    
    # Create a Polygon from the LinearRing
    poly = Polygon(ring)
    return poly

def main():
    outFileDir = os.getcwd()  # Use current directory
    p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
    p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
    n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)

    levelFacets = {}
    QTMID = {}
    levelFacets[0] = []
    QTMID[0] = []
    geojson_features = []  # Store GeoJSON features separately

    outFileName = f"octahedron.geojson"
    outFile = os.path.join(outFileDir, outFileName)

    initial_facets = [
        [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
        [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
        [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
        [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
        [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
        [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
        [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
        [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
    ]

    for i, facet in enumerate(initial_facets):
        QTMID[0].append(str(i + 1))
        geojson_features.append({
            "type": "Feature",
            "geometry": constructGeometry(facet).__geo_interface__,
            "properties": {"qtm": QTMID[0][i]}
        })
        levelFacets[0].append(facet)
        

    with open(outFile, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f)

    print(f"Octahedron saved as {outFile}")

if __name__ == '__main__':
    main()
# Reference: https://github.com/paulojraposo/QTM/blob/master/qtmgenerator.py

# -*- coding: utf-8 -*-

#   .-.
#   /v\    L   I   N   U   X
#  // \\
# /(   )\
#  ^^-^^

# This script makes a Quarternary Triangular Mesh (QTM) to tessellate the world based
# on an octohedron, based on Geoffrey Dutton's conception:
#
# Dutton, Geoffrey H. "Planetary Modelling via Hierarchical Tessellation." In Procedings of the
# AutoCarto 9 Conference, 462–71. Baltimore, MD, 1989.
# https://pdfs.semanticscholar.org/875e/12ce948012b3eced58c0f1470dba51ef87ef.pdf
#
# This script written by Paulo Raposo (pauloj.raposo [at] outlook.com) and Randall Brown
# (ranbrown8448 [at] gmail.com). Under MIT license (see LICENSE file).
#
# Dependencies: 
#   - nvector (see https://pypi.python.org/pypi/nvector and http://www.navlab.net/nvector),
#   - OGR Python bindings (packaged with GDAL).

import os, argparse, logging, datetime, math
from shapely.geometry import Polygon, LinearRing
import json

def findCrossedMeridiansByLatitude(vert1, vert2, newLat):

    """For finding pair of meridians at which a great circle defined by two points crosses the given latitude."""
    
    # Credit to Chris Veness, https://github.com/chrisveness/geodesy.

    theta = math.radians(newLat)

    theta1 = math.radians(vert1[0])
    lamb1 = math.radians(vert1[1])
    theta2 = math.radians(vert2[0])
    lamb2 = math.radians(vert2[1])

    dlamb = lamb2 - lamb1

    x = math.sin(theta1) * math.cos(theta2) * math.cos(theta) * math.sin(dlamb)
    y = math.sin(theta1) * math.cos(theta2) * math.cos(theta) * math.cos(dlamb) - math.cos(theta1) * math.sin(theta2) * math.cos(theta)
    z = math.cos(theta1) * math.cos(theta2) * math.sin(theta) * math.sin(dlamb)

    if (z*z > x*x + y*y):
         print("Great circle doesn't reach latitude.")

    lambm = math.atan2(-y, x)
    dlambI = math.acos(z / math.sqrt(x*x + y*y))

    lambI1 = lamb1 + lambm - dlambI
    lambI2 = lamb1 + lambm + dlambI

    lon1 = (math.degrees(lambI1) + 540) % 360-180
    lon2 = (math.degrees(lambI2) + 540) % 360-180

    return lon1, lon2

def lonCheck(lon1, lon2, pointlon1, pointlon2):
    lesser, greater = sorted([pointlon1, pointlon2])
    if lon1 > lesser and lon1 < greater:
        return lon1
    else:
        return lon2

def GetMidpoint(vert1, vert2):
    midLat = (vert1[0] + vert2[0]) / 2
    midLon = (vert1[1] + vert2[1]) / 2
    return(float(midLat), float(midLon))

def constructGeometry(facet):
    """Accepting a list from this script that stores vertices, return a Shapely Polygon object."""
    if len(facet) == 5:
        # This is a triangle facet of format (vert,vert,vert,vert,orient)
        vertexTuples = facet[:4]
    if len(facet) == 6:
        # This is a rectangle facet of format (vert,vert,vert,vert,vert,northboolean)
        vertexTuples = facet[:5]

    # Create a LinearRing with the vertices
    ring = LinearRing([(vT[1], vT[0]) for vT in vertexTuples])  # sequence: lon, lat (x,y)
    
    # Create a Polygon from the LinearRing
    poly = Polygon(ring)
    return poly

def divideFacet(aFacet):
    """Will always return four facets, given one, rectangle or triangle."""

    # Important: For all facets, first vertex built is always the most south-then-west, going counter-clockwise thereafter.

    if len(aFacet) == 5:

        # This is a triangle facet.

        orient = aFacet[4]  # get the string expressing this triangle's orientation

        #       Cases, each needing subdivision:
        #                    ______           ___   ___
        #       |\      /|   \    /   /\     |  /   \  |     ^
        #       | \    / |    \  /   /  \    | /     \ |     N
        #       |__\  /__|     \/   /____\   |/       \|
        #
        #        up    up     down    up     down    down   -- orientations, as "u" or "d" in code below.


        # Find the geodetic bisectors of the three sides, store in sequence using edges defined
        # by aFacet vertex indeces: [0]&[1] , [1]&[2] , [2]&[3]
        newVerts = []

        for i in range(3):
            if aFacet[i][0] == aFacet[i+1][0] or aFacet[i][1] == aFacet[i+1][1]:
                newVerts.append(GetMidpoint(aFacet[i], aFacet[i+1]))
            else:
                newLat = (aFacet[i][0] + aFacet[i+1][0]) / 2
                newLon1, newLon2 = findCrossedMeridiansByLatitude(aFacet[i], aFacet[i + 1], newLat)

                newLon = lonCheck(newLon1, newLon2, aFacet[i][1], aFacet[i+1][1])

                newVert = (newLat, newLon)
                newVerts.append(newVert)

        if orient == "u":
            #          In the case of up facets, there will be one "top" facet
            #          and 3 "bottom" facets after subdivision; we build them in the sequence inside the triangles:
            #
            #                   2
            #                  /\         Outside the triangle, a number is the index of the vertex in aFacet,
            #                 / 1\        and a number with an asterisk is the index of the vertex in newVerts.
            #             2* /____\ 1*
            #               /\ 0  /\
            #              /2 \  /3 \
            #             /____\/____\
            #           0or3   0*     1

            newFacet0 = [newVerts[0], newVerts[1], newVerts[2], newVerts[0], "d"]
            newFacet1 = [newVerts[2], newVerts[1], aFacet[2], newVerts[2], "u"]
            newFacet2 = [aFacet[0], newVerts[0], newVerts[2], aFacet[0], "u"]
            newFacet3 = [newVerts[0], aFacet[1], newVerts[1], newVerts[0], "u"]

        if orient == "d":
            #          In the case of down facets, there will be three "top" facets
            #          and 1 "bottom" facet after subdivision; we build them in the sequence inside the triangles:
            #
            #            2_____1*_____1
            #             \ 2  /\ 3  /
            #              \  / 0\  /    Outside the triangle, a number is the index of the vertex in aFacet,
            #               \/____\/     and a number with an asterisk is the index of the vertex in newVerts.
            #              2*\ 1  /0*
            #                 \  /
            #                  \/
            #                 0or3

            newFacet0 = [newVerts[2], newVerts[0], newVerts[1], newVerts[2], "u"]
            newFacet1 = [aFacet[0], newVerts[0], newVerts[2], aFacet[0], "d"]
            newFacet2 = [newVerts[2], newVerts[1], aFacet[2], newVerts[2], "d"]
            newFacet3 = [newVerts[0], aFacet[1], newVerts[1], newVerts[0], "d"]

    if len(aFacet) == 6:

        # This is a rectangle facet.

        northBoolean = aFacet[5]  # true for north, false for south

        if northBoolean:

            # North pole rectangular facet.

            # Build new facets in the sequence inside the polygons:

            #          3..........2   <-- North Pole
            #           |        |
            #           |   1    |    Outside the polys, a number is the index of the vertex in aFacet,
            #           |        |    and a number with an asterisk is the index of the vertex in newVerts.
            #           |        |
            #         2*|--------|1*           /\
            #           |\      /|  on globe  /__\
            #           | \ 0  / |  -------> /\  /\
            #           |  \  /  |          /__\/__\
            #           | 2 \/ 3 |
            #       0or4''''''''''1
            #               0*

            newVerts = []

            for i in range(4):
                if i != 2:
                    # on iter == 1 we're going across the north pole - don't need this midpoint.

                    if aFacet[i][0] == aFacet[i+1][0] or aFacet[i][1] == aFacet[i+1][1]:
                        newVerts.append(GetMidpoint(aFacet[i], aFacet[i+1]))
                    else:
                        newLat = (aFacet[i][0] + aFacet[i+1][0])/2
                        newLon1, newLon2 = findCrossedMeridiansByLatitude(aFacet[i], aFacet[i + 1], newLat)

                        newLon = lonCheck(newLon1, newLon2, aFacet[i][1], aFacet[i+1][1])

                        newVert = (newLat, newLon)
                        newVerts.append(newVert)

            newFacet0 = [newVerts[0], newVerts[1], newVerts[2], newVerts[0], "d"]  # triangle
            newFacet1 = [newVerts[2], newVerts[1], aFacet[2], aFacet[3], newVerts[2], True]  # rectangle
            newFacet2 = [aFacet[0], newVerts[0], newVerts[2], aFacet[0], "u"]  # triangle
            newFacet3 = [newVerts[0], aFacet[1], newVerts[1], newVerts[0], "u"]  # triangle

        else:

            # South pole rectangular facet

            #               1*
            #          3..........2
            #           | 2 /\ 3 |     Outside the polys, a number is the index of the vertex in aFacet,
            #           |  /  \  |     and a number with an asterisk is the index of the vertex in newVerts.
            #           | / 0  \ |
            #           |/      \|           ________
            #         2*|--------|0*         \  /\  /
            #           |        |  on globe  \/__\/
            #           |   1    |  ------->   \  /
            #           |        |              \/
            #           |        |
            #       0or4'''''''''1   <-- South Pole

            newVerts = []

            for i in range(4):
                if i != 0:
                    # on iter == 3 we're going across the south pole - don't need this midpoint
                    if aFacet[i][0] == aFacet[i+1][0] or aFacet[i][1] == aFacet[i+1][1]:
                        newVerts.append(GetMidpoint(aFacet[i], aFacet[i+1]))
                    else:
                        newLat = (aFacet[i][0] + aFacet[i+1][0])/2
                        newLon1, newLon2 = findCrossedMeridiansByLatitude(aFacet[i], aFacet[i + 1], newLat)

                        newLon = lonCheck(newLon1, newLon2, aFacet[i][1], aFacet[i+1][1])

                        newVert = newLat, newLon
                        newVerts.append(newVert)

            newFacet0 = [newVerts[2], newVerts[0], newVerts[1], newVerts[2], "u"]  # triangle
            newFacet1 = [aFacet[0], aFacet[1], newVerts[0], newVerts[2], aFacet[0], False]  # rectangle
            newFacet2 = [newVerts[2], newVerts[1], aFacet[3], newVerts[2], "d"]  # triangle
            newFacet3 = [newVerts[1], newVerts[0], aFacet[2], newVerts[1], "d"]  # triangle

    # In all cases, return the four facets made in a list
    return [newFacet0, newFacet1, newFacet2, newFacet3]


def printandlog(msg):

    """Given a string, this will both log it and print it to the console."""
    print(msg)
    logging.info(msg)


def main():
    parser = argparse.ArgumentParser(description='Builds a Dutton QTM and outputs it as GeoJSON files in WGS84 coordinates.')
    parser.add_argument('-r', '--resolution', help='Resolution level to generate.', required=True, type=int)
    args = parser.parse_args()

    resolution = args.resolution
    outFileDir = os.getcwd()  # Use current directory

    logging.basicConfig(filename=os.path.join(outFileDir, "qtm_creation_log.txt"), level=logging.DEBUG)
    startTime = datetime.datetime.now()
    printandlog(f"Starting at {startTime}")

    # Build the vertices of the initial 8 octohedron facets, as rectangles. Tuples: (Lat,Lon)
    #
    # N Pole  ---> (90,-180) ------- (90,-90)  ------- (90,0)  -------- (90, 90)  ------ (90, 180)
    #                  |                 |                |                |                 |
    #                  |                 |          Prime Meridian         |                 |
    #                  |                 |                |                |                 |
    # Equator ---> (0, -180)  ------ (0, -90)  ------- (0,0)  ---------  (0, 90)  ------ (0, 180)
    #                  |                 |                |                |                 |
    #                  |                 |                |                |                 |
    #                  |                 |                |                |                 |
    # S Pole  ---> (-90, -180) ----- (-90, -90) ------ (-90,0) -------- (-90, 90) ------ (-90, 180)
    
    p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
    p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
    n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)

    levelFacets = {}
    QTMID = {}

    for lvl in range(resolution):
        levelFacets[lvl] = []
        QTMID[lvl] = []
        geojson_features = []  # Store GeoJSON features separately

        outFileName = f"qtm_level_{lvl + 1}.geojson"
        outFile = os.path.join(outFileDir, outFileName)

        if lvl == 0:
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
                    "properties": {"QTMID": QTMID[0][i]}
                })
                levelFacets[0].append(facet)
        else:
            for i, pf in enumerate(levelFacets[lvl - 1]):
                subdivided_facets = divideFacet(pf)
                for j, subfacet in enumerate(subdivided_facets):
                    new_id = QTMID[lvl - 1][i] + str(j)
                    QTMID[lvl].append(new_id)
                    geojson_features.append({
                        "type": "Feature",
                        "geometry": constructGeometry(subfacet).__geo_interface__,
                        "properties": {"QTMID": new_id}
                    })
                    levelFacets[lvl].append(subfacet)

        with open(outFile, 'w') as f:
            json.dump({"type": "FeatureCollection", "features": geojson_features}, f)

        printandlog(f"GeoJSON for level {lvl + 1} saved as {outFile}")

    printandlog(f"Finished at {datetime.datetime.now()}")

if __name__ == '__main__':
    main()
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

################ 
# Modified by Vgrid
################

import math
from shapely.geometry import Point, Polygon, LinearRing

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

def qtm_id_to_facet(qtm_id):    
    # Base octahedral face definitions (these define the initial 8 triangular faces of QTM)
    p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
    p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
    n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)

    initial_facets = {
            '1': [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
            '2': [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
            '3': [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
            '4': [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
            '5': [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
            '6': [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
            '7': [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
            '8': [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
        }
    
    base_facet = initial_facets.get(qtm_id[0])
    if base_facet is None:
        raise ValueError("Invalid QTM ID: Base facet must be 1-8")

    # Recursively follow the QTM ID to refine the facet
    facet = base_facet
    for level in range(1, len(qtm_id)):
        facet = divideFacet(facet)[int(qtm_id[level])]

    return facet

def latlon_to_qtm_id(lat, lon, resolution):   
    # Base octahedral face definitions
    p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
    p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
    n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)

    # Initial 8 facets
    initial_facets = {
        '1': [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
        '2': [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
        '3': [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
        '4': [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
        '5': [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
        '6': [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
        '7': [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
        '8': [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
    }

    # Find the initial facet containing (lat, lon)
    for facet_id, facet in initial_facets.items():
        facet_geom = constructGeometry(facet)  
        if facet_geom.contains(Point(lon, lat)):
            qtm_id = facet_id
            current_facet = facet
            break
    else:
        raise ValueError("Point is outside the valid range")

    # Refine the facet through subdivisions up to the given resolution
    for _ in range(1, resolution):
        subfacets = divideFacet(current_facet)
        for i, subfacet in enumerate(subfacets):
            subfacet_geom = constructGeometry(subfacet)  
            if subfacet_geom.contains(Point(lon, lat)):
                qtm_id += str(i)
                current_facet = subfacet
                break

    return qtm_id

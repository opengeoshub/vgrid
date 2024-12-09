"""
GeoTransform functions relating to EASE-DGGS.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
"""
# © Regents of the University of Minnesota. All rights reserved.
# This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.

import math

import numpy as np
import rasterio
from geopandas import GeoSeries
from itertools import chain, product
from pyproj import Transformer
from shapely.geometry import Point

from vgrid.utils.easedggs.constants import levels_specs, grid_spec, ease_crs, geo_crs

from vgrid.utils.easedggs.dggs.grid_addressing import coords_ease_to_coords_grid, _grid_xy_to_grid_id
from vgrid.utils.easedggs.dggs.utils import shift_range_grid_multiple, pairwise_circle, flatten, \
    add_nodes, epsilon_check

######
#
# Helper functions
#
######
# define the pyproj Transforms necessary. need both a
#   from WGS84 -> EASE v2 (w2e) and
#   from EASE v2 -> WGS84 (e2w)
w2e = Transformer.from_crs(geo_crs, ease_crs, always_xy=True).transform
e2w = Transformer.from_crs(ease_crs, geo_crs, always_xy=True).transform


def grid_id_to_corner_coord(gid, level, shift=False):
    """Convert EASE-DGGS ID into an EASE v2 corner coordinate.

    Parameters
    ----------
    gid : str
        Grid ID of the cell.
    corner : str
        The type of corner to generate coordinates for
        (upper_left | upper_right | lower_right | lower_left)
    Returns
    ----------
    corner_x, corner_y : np.array
        EASE grid coordinate (x, y) of the cell' corner. Resolution determined using grid ID.
    """
    # r_digit = 6
    # decompse the EASE-DGGS ID into its level parts. Then remove the first compoenent
    #   which, is simply represent the Level of the cell.
    id_elements = gid.split('.')
    # level = int(id_elements[0][-1])
    id_elements.pop(0)

    # the gridID contains the index values for number of rows|columns, from a 0,0 orgin
    #   This is the 'global' position within the larger grid
    #
    #  Remember format is Lx.RRRCCC.RC.RC ...; y=R, x=C
    #
    #   Those R, C components indicate the location/number of grid cells that
    #   correspond with the location in ease grid space
    #
    #   This can be explicitly exploited such that:
    #   ([RRR, R, ... R] * [level[0]['x_length], ... level[6]['x_length]]).sum()
    #   ([CCC, C, ... C] * [level[0]['y_length],... level[6]['y_length]]).sum()
    #   will give the coordinates to the edge of the pixel. For example:
    #
    #   [L6.202482.13.21.00.00.00.00] would decompose into
    #     x = [482, 3, 2, 0, 0, 0, 0]
    #     y = [202, 1, 1, 0, 0, 0, 0]
    #
    #   To determining the coordinate of the left, top edges at Level 2
    #     x = 482 * level[0]['x_length] + 3 * [level[1]['x_length] + 2 * [level[`2`]['x_length]
    #     y = 202 * [level[0]['x_length] + 1 * [level[1]['x_length] + 1 * [level[2]['x_length]
    #
    #   and to determin if it is an edge at level 2, all the remaining Level components must be
    #       zero: [ 0, 0, 0, 0]] summing them all ensure that they are zero
    #
    #   if we asked for the Level 1 coordinate, it would be different. x|y coord is calculated
    #   using the L0, L1 compoents:
    #   x = 482 * level[0]['x_length] + 3 * [level[1]['x_length] + 2
    #
    #   but then we need to check if it is on an L1 edge. this requires all other compontes to be
    #       zeros:
    #       [2, 0, 0, 0, 0]
    #   since it is not a Level 1 edge, if we needed a right, bottom value, we'd have to shfit
    #       the last L1 index up by one

    # beyond level 0, the row | column index is grid location within the L0 pixes
    #    This is the relative position within the gloabl position.
    # spliting these into indexes, which requires converstion from str to ints
    level_x_i = np.array([int(id[1]) if len(id) != 6 else int(id[3:]) \
                          for id in id_elements], dtype=int)
    level_y_i = np.array([int(id[0]) if len(id) != 6 else int(id[0:3]) \
                          for id in id_elements], dtype=int)

    # remainders of levels are useful for figuring out if on edge on not
    #   depending on supplied level, we can drop unnecesary levels for
    #   determining the coodinates, but only need the levels < specified
    #   level for determing if it is an edge pixel or not
    if int(level) < 6:
        level_x_r = level_x_i[int(level) + 1:]
        level_y_r = level_y_i[int(level) + 1:]
    else:
        level_x_r = level_x_i[int(level):]
        level_y_r = level_y_i[int(level):]

    level_x_i = level_x_i[:int(level) + 1]
    level_y_i = level_y_i[:int(level) + 1]

    if (level_x_r.sum() != 0) and shift:
        # print(f'level_x_r: {level_x_r}; sum: {level_x_r.sum()}')
        # print('shifting x')
        level_x_i[-1] += 1

    if (level_y_r.sum() != 0) and shift:
        # print(f'level_y_r: {level_y_r}; sum: {level_y_r.sum()}')
        # print('shifting y')
        level_y_i[-1] += 1

    # print (f'level_x_i: {level_x_i}; level_y_i: {level_y_i}')
    corner_x = np.array([level_x_i[i] * levels_specs[i]['x_length'] \
                         for i in range(level_x_i.shape[0])]).sum()

    corner_y = np.array([level_y_i[i] * levels_specs[i]['y_length'] \
                         for i in range(level_y_i.shape[0])]).sum()

    return corner_x, corner_y


######
#
# main functions
#
######


def easedggs_grid_bounds(bounds, source_crs, level):
    """Determine EASE v2 coords corresponding to bounds of EASE-DGGS

    Parameters
    ----------
    bounds : tuple
        Tuple with the (min_x/left, min_y/bottom, max_x/right, max_y/top) of the soure raster.
    level : int
        EASE-DGGS Level of the
    source_crs : int
        EPSG code of coordinate reference systems for the source raster.

    Returns
    ----------
    easedggs_grid_bounds : tuple of floats
        Bounding box coordidinates with valid EASE-DGGS coordinates.
    """
    # first, convert the bounds coordinates into corner coordinates. these are tuples.
    #   next, convert the corners to 'line_segements'. these are the coordinate pairs
    #   connecting the corners
    #
    # Those line segmensts have 'duplicate' entries, since each line segment ends where next starts
    #   once we've delt with non-ease sourc_src below, we'll drop the final element of
    #   of line_segments, to reduce duplicates
    lower_left, upper_left, lower_right, upper_right = \
        chain(product([bounds[0], bounds[2]], [bounds[1], bounds[3]]))
    # print(f'conrner passed in: {corners}')
    line_segments = list(pairwise_circle([upper_left, upper_right, lower_right, lower_left]))

    # if not in EASE v2, need to build a cood tranformer to get coords into EASE
    #   also need to add nodes to the each line segment. once that is done, we use
    #   the transformer to convert all individual nodes to EASE, while maintaining the
    #   original order.
    #   then use
    if source_crs != ease_crs:
        line_segments = [add_nodes(ls[0], ls[1]) for ls in line_segments]
        tranform_coords = Transformer.from_crs(source_crs, ease_crs, always_xy=True).transform
        line_segments = [[tranform_coords(c[0], c[1]) for c in seg] for seg in line_segments]

    # may not need to drop the last elements in each list (recall, the last element of
    #   one list, is the start of the next) since we'll just be using the nodes as Points
    #   in a geoseries. but if we wanted to make Polygon, the duplciate values pose
    #   problems
    line_segments = flatten([seg[:-1] for seg in line_segments])
    bound_box = GeoSeries([Point(c[0], c[1]) for c in line_segments],
                          crs=ease_crs)
    bound_box = bound_box.total_bounds

    # now that we have ensured we have accounted for bounds not necessarily
    #   corresponding to vertices in EASE coordinates, we'll create the
    #   bounding box as vertices
    bound_box = GeoSeries([Point(bound_box[0], bound_box[3]),  # upper_left
                           Point(bound_box[2], bound_box[1])],  # lower_right
                          crs=ease_crs)

    bbox_grid_xy = coords_ease_to_coords_grid(bound_box)
    bbox_gridID = bbox_grid_xy.apply(lambda x: _grid_xy_to_grid_id(x, 6))

    # if we use upper_right, and lower left coords only, we can apply
    #   shift or not, to both points.
    upper_left = grid_id_to_corner_coord(bbox_gridID.iloc[0], level=level)
    lower_right = grid_id_to_corner_coord(bbox_gridID.iloc[1], level=level, shift=True)

    ease_xs = [upper_left[0], lower_right[0]]
    ease_ys = [lower_right[1], upper_left[1]]

    ease_xs = np.array(
        [shift_range_grid_multiple(ease_xs[i], axis='x') \
         for i in range(len(ease_xs))]
    )
    ease_ys = np.array(
        [shift_range_grid_multiple(ease_ys[i], axis='y') \
         for i in range(len(ease_ys))]
    )

    return [ease_xs[0], ease_ys[0], ease_xs[1], ease_ys[1]]  # , decimals=8)


def easedggs_align_check(in_raster, level):
    """Check if a raster is aligned to the EASE-DGGS.

    Takes a path to a raster and performs 6 tests of alignment with the
    EASE-DGGS. The first 4 tests ensure the grid is defined the same as
    the EASE-DGGS. If a raster is the same projection (EPSG:6933), is
    rectilinear (aligned to the axes), has square cells, and is one of
    the 7 EASE-DGGS resolutions, it is able to be aligned to the EASE-DGGS.
    The final 2 tests check that if you were to extrapolate your
    raster out to a global extent, that the upper left corner of the
    EASE-DGGS falls on the upper left corner of one of the input
    raster's extraploted cells.


    Parameters
    ------------
    in_raster: str, file object or pathlib.Path object
        Filename and path of input raster.

    level: int
        Integer from 0 to 6 corresponding to the EASE-DGGS levels
        (i.e., resolutions).

    Returns
    -------
    results: Dictionary
        Dictionary of test name and result of test (pass/fail).

    """
    if level not in [0, 1, 2, 3, 4, 5, 6]:
        raise Exception("Invalid grid level. Options are: 0, 1, 2, 3, 4, 5, 6")

    # Calculate the geotransform of the EASE-DGGS
    easedggs_t = rasterio.transform.from_origin(
        grid_spec['ease']['min_x'],
        grid_spec['ease']['max_y'],
        levels_specs[level]['x_length'],
        levels_specs[level]['x_length']
    )

    with rasterio.open(in_raster, 'r') as src:
        src_t = src.transform
        test_count = 6

        results = {
            'Projection is EASE-Grid 2.0 (EPSG:6933):': 'Pass',
            'Dataset is rectilinear:': 'Pass',
            'X and Y size are equal:': 'Pass',
            'Cell size matches a EASE-DGGS resolution:': 'Pass',
            'The X dimension of the cell corners matches a EASE-DGGS corner:': 'Pass',
            'The Y dimension of the cell corners matches a EASE-DGGS corner:': 'Pass',
            'Tests passed out of 6:': test_count
        }

        if src.crs != 'EPSG: 6933':
            results['Projection is EASE-Grid 2.0 (EPSG:6933):'] = str('Fail. Supplied CRS is ' + str(src.crs))
            test_count -= 1

        if not src_t.is_rectilinear:
            results['Dataset is rectilinear:'] = 'Fail'
            test_count -= 1

        if not epsilon_check(abs(src_t.a), abs(src_t.e)):
            results['X and Y size are equal:'] = 'Fail'
            test_count -= 1

        if not epsilon_check(easedggs_t.a, src_t.a):
            results['Cell size matches a EASE-DGGS resolution:'] = 'Fail'
            test_count -= 1

        # Find the input raster coordinates of the hypothetical input raster cell that overlaps the GEMS origin
        xr, yr = ~src_t * (easedggs_t.xoff, easedggs_t.yoff)

        # Move to the coordinates of the upper left corner of the above cell
        xrnd = math.floor(xr + 0.5) if xr > 0.0 else math.ceil(xr - 0.5)
        if epsilon_check(xrnd, xr):  # Test if on an edge, else go left
            xr = float(xrnd)
        else:
            xr = float(math.floor(xr))

        yrnd = math.floor(yr + 0.5) if yr > 0.0 else math.ceil(yr - 0.5)
        if epsilon_check(yrnd, yr):  # Test if on an edge, else go up
            yr = float(yrnd)
        else:
            yr = float(math.floor(yr))

        # Convert from input raster coordinates to world (ease) coordinates
        xw, yw = src_t * (xr, yr)

        # Compare x and y values calculated in previous step to gems origin coordinates
        if not epsilon_check(xw, easedggs_t.xoff):
            results['The X dimension of the cell corners matches a EASE-DGGS corner:'] = \
                str('Fail. Difference = ' + str(round(abs(xw - easedggs_t.xoff), 5)))
            test_count -= 1

        if not epsilon_check(yw, easedggs_t.yoff):
            results['The Y dimension of the cell corners matches a EASE-DGGS corner:'] = \
                str('Fail. Difference = ' + str(round(abs(yw - easedggs_t.yoff), 5)))
            test_count -= 1

        else:
            pass
        results['Tests passed out of 6:'] = str(test_count)

        # for item in results:
        #     print(item, results[item])

        return results

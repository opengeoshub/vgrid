'''
This module deals with indexing EASE-DGGS coordinates.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''
# import re

from itertools import product

import numpy as np

# from geopandas import GeoSeries, GeoDataFrame
import geopandas as gpd
from shapely import wkt
from shapely.geometry import Point#, Polygon

from vgrid.utils.easedggs.constants import grid_spec, levels_specs, ease_crs, geo_crs, cell_scale_factors

from vgrid.utils.easedggs.dggs.utils import pairwise_circle, flatten
from vgrid.utils.easedggs.dggs.utils import format_response, gen_point_grid, get_polygon_corners

from vgrid.utils.easedggs.dggs.checks import check_level, validate_coords_lon_lat, validate_grid_ids, \
    check_coords_range
from vgrid.utils.easedggs.dggs.transforms import coords_lon_lat_to_coords_ease, coords_ease_to_coords_grid,\
    grid_xy_coord_to_ease_coord
from vgrid.utils.easedggs.logConfig import logger

######
#
# Main indexing functions
#
######

def geos_to_grid_ids(coords_lon_lat, level=0, levels_specs=levels_specs, source_crs = 4326, target_crs=ease_crs):
    '''
    Return the EASE-DGGS ID for the cell correpsonding with lon, lat pair.

    Parameters
    ----------
    coords_lon_lat : list
       List of longitude,latitude cooridnate pairs (lon, lat).
    level _ int
        The specific resolution in the heiararch for the corresponding grid cell. Default is 0, the coarsest resolution.
    levels_specs : dictonary
        The dictionary with paratmer and config options. Default is 'levels_specs'

    Returns
    -------
    grid_id : dict
        The grid IDs of the cells for the cooridnates, at the specified resolution that correspond with (lon, lat)
    '''

    if not check_level(level):
        success = False
        data = ['The specified level is invalid.']

        return format_response(data, success)

    success, data = validate_coords_lon_lat(coords_lon_lat)
    if not success:
        return format_response(data, success)

    coords_ease = coords_lon_lat_to_coords_ease(coords_lon_lat)

    coords_grid = coords_ease_to_coords_grid(coords_ease)

    grid_ids = coords_grid.apply(lambda coord: _grid_xy_to_grid_id(coord, level = level))

    return format_response(grid_ids.tolist(), success)

def _gid_to_coord_ease(gid, cell_scale_factors  = cell_scale_factors, centroid_offset = 0.5):
    '''
    Convert EASE-DGGS cell ID to ease_x, ease_y cooridnates.

    Parameters
    ----------
    grid_id : String
       EASE-DGGS cell ID.

    cell_scale_factors : numpy array
       Array with the level scaling factors .

    Returns
    -------
    ease_x, ease_y : tuple
        ease_x, ease_y cooridnates of the specified cell.
    '''
    if not isinstance(gid, str):
        return False
    # centroid are half of the grid cell x,y
    # centroid_offset = 0.5

    id_elements = gid.split('.')

    # level is used in indexing the  pixel_scale_factors array
    #   add one to get the indexing right
    index_max = int(id_elements[0][1]) + 1

    id_elements.pop(0)

    # Level 0 has 3 digits for rows, columns. all other levels have 1 each
    #   deal with Level 0, then remove it from elements. then append the list
    y_row_elements = [float(id_elements[0][0:3])]
    x_col_elements = [float(id_elements[0][3:])]

    id_elements.pop(0)

    for id in id_elements:
        y_row_elements.append(float(id[0]))
        x_col_elements.append(float(id[1]))

    y_row_elements = np.array(y_row_elements)
    x_col_elements = np.array(x_col_elements)

    # remove unused levels from the arrays; makes book keeing easier
    cell_scale_factors = cell_scale_factors[0:index_max]

    y_row_scaled = y_row_elements * cell_scale_factors
    x_col_scaled = x_col_elements * cell_scale_factors

    grid_offset = (cell_scale_factors[-1] * centroid_offset)

    x_grid = np.cumsum(x_col_scaled)[-1] + grid_offset
    y_grid = np.cumsum(y_row_scaled)[-1] + grid_offset
#     print('x_grid: {}; y_grid: {}'.format(x_grid, y_grid))
    logger.debug('x_grid: {}; y_grid: {}'.format(x_grid, y_grid))

    ease_x, ease_y = grid_xy_coord_to_ease_coord(x_grid = x_grid, y_grid = y_grid )

    return Point(ease_x, ease_y)#, grid_offset

def grid_ids_to_ease(grid_ids, cell_scale_factors  = cell_scale_factors, target_crs = ease_crs):
    '''
    Convert a list of EASE-DGGS IDs to EASE Grid v2 cooridnates: Point(ease_x, ease_y).

    Parameters
    ----------
    grid_ids : List
       List of EASE-DGGS cell IDs to convert to Point(ease_x, ease_y).

    cell_scale_factors : numpy array
       Array with the level scaling factors for each level.

    Returns
    -------
    coords_ease : GeoSeries
        GeoSeies of coordinates (ease_x, ease_y) for corresponding grid IDs.
    '''

    if not isinstance(grid_ids, list):
        return False

    coords_ease = [_gid_to_coord_ease(gid) for gid in grid_ids]

    coords_ease = gpd.GeoSeries(coords_ease, crs = target_crs)

    return coords_ease

def grid_ids_to_geos(grid_ids, cell_scale_factors  = cell_scale_factors, source_crs=ease_crs, target_crs=geo_crs):
    '''
    Convert EASE-DGGS cell ID to ease_x, ease_y cooridnates.

    Parameters
    ----------
    grid_ids : List
       List of EASE-DGGS cell IDs to convert to (lon, lat).

    cell_scale_factors : numpy array
       Array with the level scaling factors for each level.

    Returns
    -------
    coords_lon_lat : dict
        Geographic coordinates (lon, lat) corresponding withe grid ID.
    '''

    success, data = validate_grid_ids(grid_ids)
    if not success:
        return format_response(data, success)

    coords_ease = grid_ids_to_ease(grid_ids)

    coords_lon_lat = coords_ease.to_crs(target_crs)
    coords_lon_lat = coords_lon_lat.apply(lambda coord: (coord.x, coord.y))
    data = coords_lon_lat.to_list()

    return format_response(data, success)

def _grid_xy_to_grid_id(grid_xy, level=0):
    '''
    Convert a EASE-DGGS coordinate (x, y) into corresponding Grid ID for specified level;

    Parameters
    ----------
    grid_xy : Shapley point(x,y)
        EASE-DGGS cooridnate point

    level : str, optional
        EASE-DGGS level/resolution. Defaults to 0.

    Returns
    ----------
    grid_id : str
        The EASE-DGGS id corrsponding to the supplied EASE-DGGS coordinate (x, y).
    '''
    if not isinstance(grid_xy, Point):
        return False

    # small number problem. This is section comes from a discussion
    #   with David Porter & other GEMS devs. Basically, it is  way
    #   of arriving at a calculated way of determing where to round
    #
    # the calculation below effectively determins that 8 decimal places
    #   is a reasonable place to round. The assumption is the EASE grid
    #   values will be anywhere between min/max, which is the space for values
    #   with the precision we're working with from pyproj etc, the 10e15 term
    #   is used to get a sense of what a small number in EASE v2 is
    #
    # in pract, rounding to 8 decimals in the grid_xy space (basically, index space)
    #   means we can detect differences in coordinates on the order of milimeters.
    #   so far, simmple testing has meant that this should be good enough enough for
    #   determining multiples in grid_xy space using np.divmod
    #
    r_digit = 6 # 6 works; 8 causes problems with y dims
    x = np.around(grid_xy.x, decimals=r_digit)
    y = np.around(grid_xy.y, decimals=r_digit)
    np.set_printoptions(precision=6)

    # small rounding differences around 0 were causing problems.
    #   in particualr, along the left hand edge of the grid,
    #   some of the x values were negative, and slightly less
    #   than 0. because np.divmod (used for dividing the cell)
    #   uses a floor operation, those numbers that were < 0, the floor
    #   results in them being -1,not 0. This can be corrected with rounding
    #   or other solutions, but since the grid is only defined for positive x,y
    #   values, the simplist solution is for
    if x < 0.0:
        x = 0.

    if y < 0.0:
        y = 0.

    grid_cell_id = ['L{}'.format(level)]

    # for every level between 0 and supplied level, inclusive
    #   want to figure our the row/column index. see documentation
    #   on the actual scheme
    for lv in range(0, level + 1, 1):
#         lv = lv
        lv = int(lv)

        # because coordinate were transformed above, nothing needed here.
        #   both x,y integer-parts represent the major Level 0 grid cell
        #       in grid cooridnate space
        #
        #   the float-part indicate the sub-location in that pixel. later, the
        #      the mod will be multipled by the refine_ratio, which will give the
        #      location in the finer mesh grids
        x_div, x_mod = np.divmod(x, 1)
        y_div, y_mod = np.divmod(y, 1)

        if lv == 0:
            # the first Level has a 3 digit space for row, columns
            #    RRRCCC format
            rc = '{:03d}{:03d}'.format(int(y_div),int(x_div))

        else:
            # everything Level after L0 uses a single integer for row, column index
            #   RC format
            rc = '{:d}{:d}'.format(int(y_div), int(x_div))

        grid_cell_id.append('{}'.format(str(rc)))

        # each time after first, we need to multiple the mod by refine_ratio value
        x = np.around(x_mod * levels_specs[lv]['refine_ratio'], decimals=r_digit)
        y = np.around(y_mod * levels_specs[lv]['refine_ratio'], decimals=r_digit)

    seperator = '.'

    return seperator.join(grid_cell_id)

def ease_polygon_to_grid_ids(polygon_ease, level=0, source_crs = ease_crs,  levels_specs = levels_specs, wkt_geom = True):
    '''
    Identify all grid cell IDs that correspond with supplied polygon

    Parameters
    ----------
    polygon : WKT
       WKT the polygon to convert to grid cell IDs.

    level _ int
        The grid level of constituent cell IDs to return

    Returns
    -------
    Grid IDs: dict
        Grid cell IDs for all constituent cells at the specified level.
    '''
    if not check_level(level):
        success = False
        data = ['The specified level is invalid.']

        return success, data

    if not wkt_geom:
        success = False
        data = ['Expected the polygon in Well Known Text (WKT) format.']

        return format_response(data, success)
    # convert to EASE Grid
    if wkt_geom:
        polygon_ease = wkt.loads(polygon_ease)

    if not ((polygon_ease.geom_type == 'MultiPolygon') or
        (polygon_ease.geom_type == 'Polygon')):

        success = False
        data = ['The input geometry should be a Polygon or MultiPolygon.']

        return format_response(data, success)

    if source_crs != 6933:
        success = False
        data = ['The expected source crs is 6933.']

        return format_response(data, success)

    success = True

    # having issues where corner coords are not quite behaving
    #   as might be expected. cell L0.202482 is causing difficulties
    #   when trying to map the locations of conrers, I am going to add
    #   small buffer to the polygon, to try and fix the problem
    #   will add 1/2 the smallest grid cell size. cap_type = 3 is
    #   square buffer
    fin_lev = list(levels_specs)[-1]
    buff_dist = levels_specs[fin_lev]['x_length'] * 0.5
    polygon_ease_in = polygon_ease
    polygon_ease = polygon_ease.buffer(buff_dist, cap_style=3)

    # poly = poly.to_crs(target_crs)

    # get corner coordiantes in EASE
    corner_coords_ease = get_polygon_corners(polygon_ease)

    corner_coords_ease = gpd.GeoSeries(corner_coords_ease, crs = source_crs)

    corner_coords_grid = coords_ease_to_coords_grid(corner_coords_ease)
    corner_grid_ids = corner_coords_grid.apply(lambda coord: _grid_xy_to_grid_id(coord, level = level))

    corner_centroids = grid_ids_to_ease(corner_grid_ids.to_list())

    # get x, y deltas for grid points; these are the distances between ceintroids
    delta_x = levels_specs[level]['x_length']
    delta_y = levels_specs[level]['y_length']

    point_grid = gen_point_grid(upper_left = corner_centroids[0],
                                upper_right = corner_centroids[1],
                                lower_right = corner_centroids[2],
                                lower_left = corner_centroids[3],
                                delta_x = delta_x , delta_y = delta_y)


    # ok, below presumes dtaframes, but the shapes input is not
    # RETHING df use, and just use the shapely options?
    polygon_ease_in = gpd.GeoDataFrame({'name':['test_shape'], 'geometry': polygon_ease_in}, crs=source_crs)


    #points_in_poly = gpd.sjoin(point_grid, polygon_ease_in, how='left', op='within')
    points_in_poly = gpd.sjoin(point_grid, polygon_ease_in, how='left', predicate='within')
    points_in_poly = points_in_poly.dropna().reset_index()

    # once the sjoin is done, there is not a lot of point carting around the df framework.
    #   may not make much difference on the ultimate performance, but no need to cart it around
    points_in_poly = points_in_poly.geometry

    coords_grid = coords_ease_to_coords_grid(points_in_poly)
    grid_ids = coords_grid.apply(lambda coord: _grid_xy_to_grid_id(coord, level = level))

    return format_response(grid_ids.to_list(), success)

def geo_polygon_to_grid_ids(polygon_lon_lat, level=0, source_crs = geo_crs, target_crs = ease_crs, levels_specs = levels_specs, return_centroids = True, wkt_geom=True):
    '''
    Identify all grid cell IDs that correspond with the supplied polygon (lon, lat).

    Parameters
    ----------
    polygon : WKT
       WKT the polygon to convert to grid cell IDs.

    level _ int
        The grid level of constituent cell IDs to return

    source_crs : int
        The source EPSG code for the polygon. Default is 4326 (lon, lat).

    target_crs : int
        The targest EPSG code for the polygon. Default is 6933 (EASE Grid v2)

    Returns
    -------
    Grid IDs : dict
        Grid cell IDs for all constituent cells at the specified level.
    '''
# convert to EASE Grid
    if not wkt_geom:
        success = False
        data = ['Expected the polygon in Well Known Text (WKT) format.']

        return format_response(data, success)

    if wkt_geom:
       polygon_lon_lat = wkt.loads(polygon_lon_lat)

    if source_crs != 4326:
        success = False
        data = ['The expected source crs is 4326.']

        return format_response(data, success)

    corners_lon_lat = get_polygon_corners(polygon_lon_lat)
    corners_lon_lat = [(coord.x, coord.y) for coord in corners_lon_lat]

    if not check_coords_range(corners_lon_lat):
        success = False

        data = f"""Lon range is {grid_spec['geo']['min_x']} : {grid_spec['geo']['max_x']} ; lat range is {grid_spec['geo']['max_y']} : {grid_spec['geo']['min_y']}"""
        data =[data]

    polygon_lon_lat = gpd.GeoSeries(polygon_lon_lat, crs = source_crs)

    polygon_ease = polygon_lon_lat.to_crs(target_crs)
    polygon_ease = polygon_ease.geometry.values[0].wkt

    response = ease_polygon_to_grid_ids(polygon_ease, level = level)

    return response

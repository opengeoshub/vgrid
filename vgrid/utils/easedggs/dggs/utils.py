'''
Module for formatting responses for API.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''
from itertools import chain, product, tee

import numpy as np
import geopandas as gpd

from shapely import wkt
from shapely.geometry import Polygon, Point

from vgrid.utils.easedggs.constants import levels_specs, ease_crs, grid_spec
from vgrid.utils.easedggs.logConfig import logger

def format_response(data, success):
    '''
    Format the library repsonse for return to API
    Args:
         data: list
            Contains data, or the error message
        success :(boolean)
            Indicatication of function success or failure.
    Returns:
        response : dict
            The formatted response to return to the API.
    '''
    if success:
        response =  {'success' : success,
                    'result' : {'data' :data}}

    else:
        response =  {'success' : success,
                    'result' : {'error_message':data}}

    return response

def boundbox_to_poly(left, bottom, right, top):
    '''
    convert bounding extents to polygon
    Parameters
    ----------
    left : float
        Left edge of bounding box.
    bottom : float
        Bottom edge of bounding box.
    right : float
        Right edge of bounding box.
    top : TYPE
        Top edge of bounding box.

    Returns
    -------
    Shapely Polygon corresponding with bounding box.
    '''
    if type(left) is int:
        left = float(left)
        bottom = float(bottom)
        right = float(right)
        top = float(top)

    ul = Point(left, top)
    ur = Point(right, top)
    lr= Point(right, bottom)
    ll = Point(left, bottom)

    poly = Polygon([ul, ur, lr, ll, ul])

    return poly

def get_polygon_corners(polygon, ccw=True):
    '''
    Get the bounds of a polygon

    Parameters
    ----------
    polygon : shapely polygon
        The polygon to get the corner coordinates for.
    ccw : boolean
        Return coordinates in counter-clockwise order. Default = False

    Returns
    -------
    bounds : tuple
        Corner coodinates for the bounding box/evevlope of the polygon.
    '''

    if type(polygon) is str:
        polygon = wkt.loads(polygon)

    # .bounds method returns (minx, miny, maxx, maxy) by default.
    #    that is the same as (left, bottom, right, top).
    #  Keep that convetion, but allow for return in clockwise order.
    #     when false: (left, top, right, bottom) or: (minx, maxy, maxx, miny)
    minx, miny, maxx, maxy = polygon.bounds

    upper_left = Point(minx, maxy)
    upper_right = Point(maxx, maxy)
    lower_left = Point(minx, miny)
    lower_right = Point(maxx, miny)

    if not ccw:
        return upper_left, lower_left, lower_right, upper_right

    return upper_left, upper_right, lower_right, lower_left

def enumerate_id_elements(refine_ratio):

    if not isinstance(refine_ratio, int):
        return False

    row_ids = [str(rr) for rr in range(refine_ratio)]
    col_ids = row_ids

    sub_elements = list(chain(product(row_ids, col_ids)))
    sub_elements = [''.join(se) for se in sub_elements]

    return sub_elements

def enumerate_grid_table_rows(x_coords, y_coords, level=0, parent_id = None):
    '''
    Generate grid polygons for cells using x, y cooridnate vectors.
    Parameters
    ----------
    x_coords : float
       X coordinate vector in EASE grid coordinates.
    y_coords : float
        Y coordinate vector in EASE grid coordinates.
    level _ int
        The grid level for the corresponding geometries.
    parent_id : String
        Denotes the cell id of the parent cell. Default is None.

    Returns
    -------
    Relevant child cell characteritics: Row ID, Column ID, Grid ID, Geometry and Centroid
    '''
    geoms = []
    r_ind = []
    c_ind = []
    grid_ids = []
    centroid = []
    if level == 0:
        neighbors = []
#         children = []

    x_i = np.arange(1,len(x_coords))
    y_i = np.arange(1,len(y_coords))

    for r,c in product(np.arange(1,len(y_coords), 1), np.arange(1,len(x_coords), 1)):
        c = int(c)
        r = int(r)

        poly = boundbox_to_poly(left = x_coords[c-1],
                                bottom = y_coords[r],
                                right = x_coords[c],
                                top = y_coords[r-1])
        geoms.append(poly.wkt)
        centroid.append(poly.centroid.wkt)
        ci = c - 1
        ri = r -1

        r_ind.append(ri)
        c_ind.append(ci)

        if level == 0:
            grid_id = 'L0.{:03d}{:03d}'.format(ri, ci)

        else:
            grid_id = '{}.{:d}{:d}'.format(parent_id, ri, ci)

        grid_ids.append(grid_id)

    return r_ind, c_ind, grid_ids, geoms, centroid

def gen_point_grid(upper_left, upper_right, lower_right, lower_left, delta_x, delta_y, target_crs = ease_crs):
    '''
    Generates a grid of equally spaced points, between geograhic extent

    upper_left : shapely Point
        cooridinate for the upper left of the geographic extent

    upper_right : shapely Point
        cooridinate for the upper right of the geographic extent

    lower_right : shapely Point
        cooridinate for the lower right of the geographic extent

    lower_left : shapely Point
        cooridnate for the lower left of the geographic extent

    delta_x : float
        horizontal distance between subsequent cooridnates

    delta_y : float
        vertical distance between subsequent cooridnates

    returns:
    point_grid : pandas GeoDataFame
        grid of evenly spaced coordinates between the geographic extent
    '''
    if not (isinstance(upper_left, Point) and isinstance(upper_right, Point) and
            isinstance(lower_right, Point) and isinstance(lower_left, Point) ):
        return False

    if not (isinstance(delta_x, float) and isinstance(delta_y, float)):
        return False

    maxx = np.max([upper_left.x, upper_right.x, lower_right.x, lower_left.x])
    minx = np.min([upper_left.x, upper_right.x, lower_right.x, lower_left.x])
    # print(maxx, minx)
    logger.debug('gen_point_grid: maxx - {}; minx - {}'.format(maxx, minx))

    maxy = np.max([upper_left.y, upper_right.y, lower_right.y, lower_left.y])
    miny = np.min([upper_left.y, upper_right.y, lower_right.y, lower_left.y])

    x_steps = np.rint((maxx - minx) / delta_x).astype('int') + 1
    y_steps = np.rint((maxy - miny) / delta_y).astype('int') + 1

    x_coords = np.linspace(start = minx, stop = maxx, num = x_steps, endpoint=True)
    y_coords = np.linspace(start = maxy, stop = miny, num = y_steps, endpoint=True)

    grid = np.meshgrid(x_coords, y_coords, indexing='xy')

    x_points = np.ravel(grid[0])
    y_points = np.ravel(grid[1])

    coords = [Point(x_points[i], y_points[i]) for i in range(x_points.shape[0])]

    point_grid = gpd.GeoDataFrame({'ease_x':x_points, 'ease_y':y_points, 'geometry' : coords}, crs=target_crs)

    return point_grid

def calc_grid_coord_vectors(min_x, min_y, max_x, max_y, level=0, levels_specs = levels_specs,
                            x_ascend = True, y_ascend = False,  ):
    '''
    Calculate pixel corner coordinates for the specified level using grid bounding coordinates

    Parameters
    ----------
    min_x : float
        EASE Grid cooridinates of left edge of bounding box.
    min_y : float
        EASE Grid cooridinates of bottom edge of bounding box.
    max_x : float
        EASE Grid cooridinates of right edge of bounding box.
    max_y : TYPE
        EASE Grid cooridinates of top edge of bounding box.

    Returns
    -------
    x_coord_vector, y_coord_vector : np.array
        Numpy arrays of x, y cooridinates at specified level.
    '''
    # user in puts the number of columns/rows of the grid.
    # number of grid edges is column/rows +1
    #
    # for level 0, coordinates are generated using the number of rows/columns for the grid
    if level == 0:
        x_col = levels_specs[level]['n_col'] + 1
        y_row = levels_specs[level]['n_row'] + 1

    # for all levels > 0, the refine_ratio is used to derive the number of cells.
    #    the refine_ratio for levels is associuated with the partent in the dict
    else:
        level = level -1
        level = str(level)
        x_col = levels_specs[level]['refine_ratio'] + 1
        y_row = levels_specs[level]['refine_ratio'] + 1

    if x_ascend:
        x_coord_vector = np.linspace(start = min_x, stop = max_x, num = x_col)
    else:
        x_coord_vector = np.linspace(start = max_y, stop = min_x, num = x_col)

    if y_ascend:
        y_coord_vector = np.linspace(start = min_y, stop = max_y, num = y_row)
    else:
        y_coord_vector = np.linspace(start = max_y, stop = min_y, num = y_row)

    return x_coord_vector, y_coord_vector

######################
#
# helper function declarations
#
######################
# adapted from: https://math.stackexchange.com/questions/914823/shift-numbers-into-a-different-range
def _shift_range(val, o_start, o_end, n_start, n_end):
    '''Map value (x) from old range to new range of values.

    Parameters
    ----------
    val : numeric
        Specific numeric value to shift
    o_start : numeric
        Start value of original range
    o_end : numeric
        End value of original range
    n_start : numeric
        Start value of new range
    n_end : _type_
        End value of new range

    Returns
    -------
    numeric
        Value mapped into new range
    '''
    return (n_end - n_start) / (o_end - o_start) * (val - o_start) + n_start

def shift_range_ease(val, axis):
    '''
    Convert 1-d EASE v2 coordinate (x | y) to EASE-DGGS coordinate.

    Parameters
    ----------
    val: float
        The value to shift from EASE coords to EASE-DGGS coords.
    axis: str
        'x' or 'y'; the axis of the coords to be shifted.

        (we should be able to use an enum and pydantic or something
        to declare an acceptable set of values)

    Returns
    -------
    numeric
        Value mapped into EASE v2 bounding box range.
    '''
    if axis == 'x':
        return  _shift_range(val, grid_spec['ease'][f"min_{axis}"],
                            grid_spec['ease'][f"max_{axis}"],
                            0, levels_specs[0]['n_col'])

    if axis == 'y':
        return _shift_range(val, grid_spec['ease'][f"max_{axis}"],
                            grid_spec['ease'][f"min_{axis}"],
                            0, levels_specs[0]['n_row'])

def shift_range_grid_xy(val, axis):
    '''
    Convert 1-d EASE-DGGS coordinate (x | y) to EASE v2 grid coordinate.

    Parameters
    ----------
    val: float
        The value to shift from EASE coords to EASE-DGGS coords.
    axis: str
        'x' or 'y'; the axis of the coords to be shifted.

        (we should be able to use an enum and pydantic or something
        to declare an acceptable set of values)
    Returns
    -------
    numeric
        Value mapped into EASE-DGGS space (x|cols [0: 964]; y|rows [0: 406])
    '''
    if axis == 'x':
        return  _shift_range(val, 0, levels_specs[0]['n_col'],
                             grid_spec['ease'][f"min_{axis}"],
                             grid_spec['ease'][f"max_{axis}"])
    if axis == 'y':
        return _shift_range(val, 0, levels_specs[0]['n_row'],
                            grid_spec['ease'][f"max_{axis}"],
                            grid_spec['ease'][f"min_{axis}"])

def shift_range_grid_multiple(val, axis):
    '''
    Convert 1-d EASE-DGGS coordinate (x | y) * levels_specs[0]['x_lenth']
    to EASE v2 grid coordinate.

    Parameters
    ----------
    val: float
        The value to shift from EASE coords to EASE-DGGS coords.
    axis: str
        'x' or 'y'; the axis of the coords to be shifted.

        (we should be able to use an enum and pydantic or something to
        declare an acceptable set of values)

    Returns
    -------
    numeric
        Value mapped from (n_rows | n_cols) * cell_widths to EASE v2 bounds of
        EASE-DGGS (min|max x|y)
    '''
    if axis == 'x':
        return  _shift_range(val, 0, levels_specs[0]['n_col'] *
                            levels_specs[0]['x_length'],
                             grid_spec['ease'][f"min_{axis}"],
                             grid_spec['ease'][f"max_{axis}"])
    if axis == 'y':
        return _shift_range(val, 0, levels_specs[0]['n_row'] *
                            levels_specs[0]['y_length'],
                            grid_spec['ease'][f"max_{axis}"],
                            grid_spec['ease'][f"min_{axis}"])
# truncate function from:
#   https://stackoverflow.com/questions/42021972/truncating-decimal-digits-numpy-array-of-floats
def trunc(vals, decimals=0):
    '''Take values and truncates them at specific decimal place

    Parameters
    ----------
    values : np.array, floats
        Floats to truncate at decimal place.
    decimals : int, optional
        The decimal place to truncate the , by default 0

    Returns
    -------
    numpy.array
        Array truncated at deimal place
    '''
    return np.trunc(vals*10**decimals)/(10**decimals)

def flatten(list):
    '''Function that converts nested lists into a flat list.

    Parameters
    ----------
    list : list
        Nested list to convert to individual elements

    Returns
    -------
    flattened list
        Single list indvidual elements
    '''
    return [li for sub_list in list for li in sub_list]

# often end up doing operations going around in a 'circle', particularly
#   operations around the pairs of a bounding box. the pairwise_circle
#   is not a new idea, and hence from stack exchange
# https://stackoverflow.com/questions/36917042/pairwise-circular-python-for-loop/36927946#36927946
def pairwise_circle(iterable, reverse=False):
    '''Creates pairwise combinations from a list, with last element looping back to first.
        E.g. [UL, UR, LR, LL] -> [(UL, UR), (UR, LR), (LR, LL), (LL, UL)]
        if reverse = True
        E.g. [UL, UR, LR, LL] -> [(LL, LR), (LR, UR), (UR, UL), (UL, LL)]

    Parameters
    ----------
    iterable : _type_
        The iterable to construct pairwise combinations of.
    reverse : bool, optional
        Perform pair wise combinations in reverse order, by default False


    Returns
    -------
    iterable
        Tuples of the pairwise combinations
    '''
    if reverse:
        iterable = reversed(iterable)

    a, b = tee(iterable)
    first = next(b, None)
    return zip(a, chain(b, (first,)))

def epsilon_check(x, y, E = 1e-5) -> bool:
    '''Check that difference between 2 values is less than EPISLON.

    Parameters
    ----------
    x : float
        First value to compare
    y : float
        Second value to compare
    E : float, optional
        The value the difference between x, y should not exceed, by default EPSILON

    Returns
    -------
    bool
        True ff difference between x, y is less than, equal to EPSILON
    '''
    return abs(x-y) <= E

    # basic idea of add_node is that a line segement can be minimally defined
#   using two points. here, these are teh start, end points. The aim of this
#   function is to add addition points into the line segment at regular intervals
#   condider:
#
#   initial condition:
#   start                    end
#   O------------------------O
#
#   after add_nodes nodes =7
#   st  n1  n2  n3  n4  n5  n6  n7  ed
#   O---o---o---o---o---o---o---o---O
def add_nodes(start, end, nodes = 21):
    '''Insert nodes (coordinate pairs) between the start, end of a line segment.

    Parameters
    ----------
    start : tuple
        Starting coordinate pair of line segment to add nodes to.
    end : tuple
        Ending coordinate pair of line segment to add nodes to.
    nodes : int, optional
        Number of nodes to add to line segement, by default 21.
    # drop_end : boolean, optional
    #     Drop the final coordinate pair of the series, by default True.
    Returns
    ----------
    new line segement : list
        List with nodes added inbetween start, end.
    '''
    xs = np.linspace(start[0], end[0], nodes)
    ys = np.linspace(start[1], end[1], nodes)

    return list(zip(xs[:], ys[:]))


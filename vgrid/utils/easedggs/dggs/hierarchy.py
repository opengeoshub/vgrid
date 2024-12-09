'''
Functions related to  EASE-DGGS hierarchy.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''

import json
from itertools import product

import numpy as  np
import pandas as pd

from shapely import wkt

from vgrid.utils.easedggs.constants import levels_specs
from vgrid.utils.easedggs.dggs.checks import validate_grid_ids, check_level

from vgrid.utils.easedggs.dggs.utils import format_response, enumerate_id_elements
from vgrid.utils.easedggs.dggs.utils import enumerate_grid_table_rows
from vgrid.utils.easedggs.dggs.utils import calc_grid_coord_vectors

from vgrid.utils.easedggs.logConfig import logger

def _child_to_parent(gid, level=0):
    '''
    Determines the parent cell (coarser) of a child at the specified level.

    Parameters
    ----------
    gid : str
        Grid ID of the child cell.
    level : str
        The specified level (resolution) of the parent to return.

    Returns
    -------
    partent_id : str
        Parent ID of of the cell.
    '''
    if not isinstance(gid, str):
        return False

    id_elements = gid.split('.')
    # test that L0 not passed in
    if id_elements[0][1] == '0':
        return False

    index_max = level + 2

    parent_id =id_elements[0:index_max]
    parent_id[0] = f'L{str(level)}'

    return ('.'.join(parent_id))

def children_to_parents(children, level=0):
    '''
    Determines the parent cells (coarser) of all children at the specified level.

    Parameters
    ----------
    children : List
        Children whose parent cells you want to identify.
    level : str
        The level of the parent cells.

    Returns
    -------
    parent_ids : list
        Parent IDs of of the cells.
    '''
    if not isinstance(children, list):
        return False

    success, data = validate_grid_ids(children)

    if not success:
        return format_response(data, success)

    data = [_child_to_parent(gid, level = level) for gid in children]

    if (not any(data)):
        data = ['Invalid grid ID passed to _child_to_parent ']
        return format_response(data, False)

    return format_response(data, success)

def gen_child_geometries(parent_geometry, parent_id, child_level, wkt_geom = True):
    '''
    Generate child cell characterisitcs from parent cells.

    Parameters
    ----------
    parent_geometry : polygon
       WKT with parent polygone.
    parent_id : float
        The cell id of the parent cell.
    child_level _ int
        The grid level of the child to create
    wkt_geom : boolean
        Denotes that partent_geometry is wkt

    Returns
    -------
    row_ID, column_ID, Grid_ID, geoms, centroid    : tuple
        Relevant child cell characteritics: Row ID, Column ID, Grid ID, Geometry and Centroid
    '''
    if wkt_geom:
        parent_geometry = wkt.loads(parent_geometry)

    min_x, min_y, max_x, max_y = parent_geometry.bounds

    x_coords, y_coords = calc_grid_coord_vectors(min_x = min_x,
                                             min_y = min_y,
                                             max_x = max_x,
                                             max_y = max_y,
                                             level = child_level,
                                             x_ascend = True,
                                             y_ascend = False)

    r_ind, c_ind, grid_id, geoms, centroid = \
        enumerate_grid_table_rows(x_coords = x_coords,
                                    y_coords = y_coords,
                                    level = child_level,
                                    parent_id = parent_id)

    return(r_ind, c_ind, grid_id, geoms, centroid)

def _parent_to_children(gid, level=1):
    '''
    Determines all the children of a single parent the specified level.

    Parameters
    ----------
    gid : str
        EASE-DGGS ID of the parent cell.
    level : str
        Level (resolution) of the children cells.

    Returns
    -------
    children : list
        List of grid IDs for all the children of the parent at the supplied level.
    '''
    if not isinstance(gid, str):
        return False

    id_elements = gid.split('.')

    parent_level = len(id_elements) - 2
    parent_cell = id_elements.pop(0)
    parent_cell = '.'.join(id_elements)
    parent_cell = f'L{str(level)}.{parent_cell}'

    index_max = level

    base_id =id_elements[0:index_max]

    refine_ratio = [levels_specs[lv]['refine_ratio'] for lv in range(parent_level,index_max)]

    level_elements = [enumerate_id_elements(rr) for rr in refine_ratio]

    level_elements = list(product(*level_elements))

    level_elements = ['.'.join(le) for le in level_elements]

    children = ['{}.{}'.format(parent_cell, le) for le in level_elements]

    return(children)

def parents_to_children(grid_ids, level = 1):
    '''
    Determines all of the children cells of the parent cell for the specified level.

    Parameters
    ----------
    grid_id : list
        List of EASE-DGGS IDs for all the parent cells.
    level : str
        Level (resolution) of the children cells to return.

    Returns
    -------
    children : list
        List of EASE-DGGS IDs for children of the parent cells.
    '''
    if not isinstance(grid_ids, list):
        data = ['Input grid IDs should be list']
        return format_response(data, False)

    success, data = validate_grid_ids(grid_ids)
    if not success:
        return format_response(data, success)

    data = [_parent_to_children(gid, level=level)for gid in grid_ids]

    return format_response(data, success)


def grid_aggregate(grid_ids, grid_vals, level = 0, method = 'mean', levels_specs = levels_specs):
    '''
    Aggregate EASE-DGGS ID to coarser spatial resolution.

    Parameters
    ----------
    grid_ids : list
       List of the EASE-DGGS IDs to aggregate.

    grid_vals : list
       List of the corresponding EASE-DGGS cell values to aggregate.

    level : str
        Target resolution grid level to aggregate to.

    method : str
        The aggregation method to employ.

    Returns
    -------
    Lists with grid_ids and aggregated values lists.
    '''

    if not isinstance(grid_ids, list) and not isinstance(grid_vals, list):
        success = False
        data = 'Lists expected for grid_ids and grid_vals.'
        return format_response(data, success)

    if not check_level(level):
        success = False
        data = ['The specified level is invalid.']

        return format_response(data, success)

    if not all(isinstance(x, (int, float)) for x in grid_vals):
        success = False
        data = 'grid_vals must be ints or floats'

        return format_response(data, success)

    allow_methods = ['count', 'first', 'last', 'mean', 'median', 'min',
                    'max', 'std', 'sum', 'var', 'prod', 'mode']

    numeric_true = ['first', 'last', 'mean', 'median', 'min',
                    'max', 'std', 'sum', 'var', 'prod']

    if method not in allow_methods:
        success = False
        data = 'Invalid aggregation method supplied.'
        return format_response(data, success)

    parent_ids = children_to_parents(grid_ids, level = level)
    parent_ids = parent_ids['result']['data']

    df = pd.DataFrame({'children': grid_ids,
                        'grid_vals': grid_vals,
                        'parents' : parent_ids})

    if method in numeric_true:
        agg = df.groupby('parents').aggregate(method, numeric_only=True)
        out_ids = agg.index.to_list()
        out_vals = agg.grid_vals.to_list()
    elif method == 'mode':
        agg = df.groupby('parents').grid_vals.apply(lambda x: x.mode()[0])
        out_ids = agg.index.to_list()
        out_vals = list(agg.values)
    elif method == 'count':
        agg = df.groupby('parents').grid_vals.apply(lambda x: x.count())
        out_ids = agg.index.to_list()
        out_vals = list(agg.values)
    else:
        logger.warning('Error in heirarchy.py. Else condition encounter that shoudl not be possible.')

    return format_response({'grid_ids' : out_ids, 'values' : out_vals}, True)


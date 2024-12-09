'''
Validation functions of EASE-DGGS DGGS.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''

from vgrid.utils.easedggs.logConfig import logger #, debug_logger
import re
from vgrid.utils.easedggs.constants import grid_spec, levels_specs

# compile re in advance, to make grid_id_to_geo more flexible
# so that you don't have to have all gids of same levelo
grid_id_regex = ['L[0-6]','[0-4][0-9]{2}[0-9]{3}', '[0-3]{2}',
                '[0-2]{2}', '[0-2]{2}', '[0-9]{2}', '[0-9]{2}', '[0-9]{2}']
gid_regexs = [re.compile(r'(^{}$)'.format('\.'.join(grid_id_regex[0:i]))) for i in range(2, 9)]

######################
#
# Error checking/ validation  functions
#
######################
def check_coord_range(coord, grid_spec=grid_spec):
    '''
    Check that a single lon, lat coordinate pair is within specified ranges.

    coords : tuple
        The tuple with a coodinate pair (lon, lat) to test.
    grid_spec : dict
        The dictionary specifiying the min/max lon, lat values.

    Returns:
    result : boolean
        It the lon, lat coordinate pair within specified ranges?
    '''
    min_lon = grid_spec['geo']['min_x']
    max_lon = grid_spec['geo']['max_x']

    min_lat = grid_spec['geo']['min_y']
    max_lat = grid_spec['geo']['max_y']

    # print(coord)
    logger.debug('check_coord_range - coord : {}'.format(coord))
    if not isinstance(coord, tuple) and not isinstance(coord, list):
        return False

    # if not ((coord[0] > min_lon) and (coord[0] <= max_lon) and
            # (coord[1] >= min_lat) and (coord[1] < max_lat)):
    if not ((coord[0] >= min_lon) and (coord[0] <= max_lon) and
            (coord[1] >= min_lat) and (coord[1] <= max_lat)):
        return False
    else:
        return True

    if not all(valid):
        return False
    else:
        return True

def check_level(level, valid_levels = list(levels_specs.keys())):
    '''
    Test if input level is valid

    Parameters
    ----------
    level : int
        The level of the grid hierarchy to test.
    valid_levels : list
        List of valid levels. Default = list is keys from from levels_specs{}

    Returns
    -------
    results : boolean
        Is the specified level valid?
    '''
    if level not in valid_levels:
        return False
    else:
        return True

def check_gid_format(gid, gid_regexs = gid_regexs):
    '''
    Test if format of supplied gid is valid.

    Parameters
    ----------
    gid : string
        The GESM grid id to test.
    gid_regexs : list
        List of precompiled regexes for testing, one per grid level

    Returns
    -------
    results : boolean
        Is the grid id format valid?
    '''

    if not isinstance(gid, str):
        return False

    id_elements = gid.split('.')

    level = int(id_elements[0][-1])
    re = gid_regexs[level]

    if not re.match(gid):
        return False
    else:
        return True

def check_gid_l0_index(gid, levels_specs=levels_specs):
    '''
    Check that a single grid ID matches specification.

    Parameters
    ----------
    grid_id: str
        The grid ID to test.
    levels_specs : dict
        Dictionary containing grid specification. Default is the master leve_dict.

    Returns:
    result : boolean
        Does the supplied grid id match the spec?
    '''
    if not isinstance(gid, str):
        return False

    # get the individual id elements, don't need the Level compenent
    # id_elements = grid_id.split('.')
    l0_element = gid.split('.')[1]

    row_ind = int(l0_element[:-3])
    col_ind = int(l0_element[3:])

    row_max =levels_specs[0]['n_row'] - 1
    col_max =levels_specs[0]['n_col'] - 1
    # print(row_ind, col_ind)
    # print(row_max, col_max)

    if not ((row_ind >= 0) and (row_ind <= row_max) and
            (col_ind >= 0) and (col_ind <= col_max)):
        return False
    else:
        return True

def check_gid_num_element(grid_id):
    '''
    Check that the length of grid ID matches its level.

    Parameters
    ----------
    grid_id: str
        The grid ID to test.

    Returns:
    result : boolean
        Is the supplied grid id valid?
    '''
    if not isinstance(grid_id, str):
        return False

    # delta; used to for comparing the length of the split grid id string,
    #   matches the spcified level. The level, lenght split is:
    #       Level 0 : 2 parts
    #       Level 1 : 3 parts
    #       Level 2 : 4 parts....
    base = 2
    id_elements = grid_id.split('.')

    # print(id_elements[0])
    level = int(id_elements[0][-1])

    #if (level + base) is not len(id_elements):
    if (level + base) != len(id_elements):
        return False
    else:
        return True

def check_gid_start_with(grid_id):
    '''
    Check that a single grid ID starts with 'L'.

    Parameters
    ----------
    grid_id: str
        The grid ID to test.

    Returns:
    result : boolean
        Does the supplied grid id start with 'L'?
    '''
    if not isinstance(grid_id, str):
        return False

    if not grid_id.startswith('L'):
        return False
    else:
        return True

####
#start of plural tests
####
def check_coords_range(coords_lon_lat, grid_spec = grid_spec):
    '''
    Check that all lon, lat coordinate pair are within specified ranges.

    coords_lon_lat : list
        The list of coodinate pairs to test.

    Returns:
    result : boolean
        Are all lon, lat coordinate pairs within specified ranges?
    '''
    if not isinstance(coords_lon_lat, list):
        return False

    valid = [check_coord_range(coord) for coord in coords_lon_lat]

    if not all(valid):
        return False
    else:
        return True

def check_grid_ids_format(grid_ids, gid_regexs = gid_regexs):
    '''
    Check that all Grid IDs indicies are within specified ranges.

    grid_ids : list
        The list of grid IDs to test.

    Returns:
    result : boolean
        Do all Grid IDs have indicies within specified ranges?
    '''
    valid = [check_gid_format(gid) for gid in grid_ids]

    if not all(valid):
        return False
    else:
        return True

def check_grid_ids_l0_index(grid_ids, levels_specs=levels_specs):
    '''
    Determines if a list of grid ids have correct index reanges for L0,

    Parameters
    ----------
    grid_ids : list
        The list of grid IDs to test.

    Returns:
    result : boolean
        Do all the grid IDs have correct ranges for L0 indicies?
    '''
    if not isinstance(grid_ids, list):
        return False

    valid = [check_gid_l0_index(gid) for gid in grid_ids]
    # print(valid)
    if not all(valid):
        #prob_ids = compress(grid_ids, not(all_valid))
        return False #, prob_ids
    else:
        return True

def check_grid_ids_num_elements(grid_ids):
    '''
    Determines if a list of grid ids have correct number of elements.

    Parameters
    ----------
    grid_ids : list
        The list of grid IDs to test.

    Returns:
    result : boolean
        Do all the grid IDs have the correct number of elements?
    '''
    if not isinstance(grid_ids, list):
        return False

    valid = [check_gid_num_element(gid) for gid in grid_ids]
    # print(valid)
    if not all(valid):
        #prob_ids = compress(grid_ids, not(all_valid))
        return False #, prob_ids
    else:
        return True

def check_grid_ids_starts_with(grid_ids):
    '''
    Check that all Grid IDs in the list start with 'L'.

    grid_ids : list
        The list of grid IDs to test.

    Returns:
    result : boolean
        Are all the grid IDs the supplied grid id valid?
    '''

    valid = [check_gid_start_with(gid) for gid in grid_ids]

    if not all(valid):
        return False
    else:
        return True

def validate_coords_lon_lat(coords_lon_lat):
    '''
    Run all the steps to check that lon, lat coordinate pairs.

    coords_lon_lat : list
        The list of coorindate pairs (lon, lat) to check.

    Returns:
    boolean, data
        Did all the (lon, lat) cooridnate pairs in the list pass all the validation checks?
        data = error message
    '''
    if not check_coords_range(coords_lon_lat):
        success = False

        data = f"""Lon range is {grid_spec['geo']['min_x']} :
                {grid_spec['geo']['max_x']} ; lat range is {grid_spec['geo']['max_y']} :
                {grid_spec['geo']['min_y']}"""
        data =[data]

        return success, data

    if not isinstance(coords_lon_lat, list):
        success = False
        data = ['coords_lon_lat is not a list']

        return success, data

    return True, None

def validate_grid_ids(grid_ids):
    '''
    Run all the steps to check that grid IDs are valid.

    grid_ids : list
        The list of grid IDs to test.

    Returns:
    boolean, data
        Did all the grid IDs in the list pass all the validation checks?
        data = error message
    '''
    if not check_grid_ids_starts_with(grid_ids):
        success = False
        data = ['Grid IDs must start with \'L\'']

        return success, data

    if not check_grid_ids_num_elements(grid_ids):
        success = False
        data = ['Grid IDs contain incorrect number of elements.']

        return success, data

    if not check_grid_ids_format(grid_ids):
        success = False
        data = ['Grid IDs contain improperly formatted IDs']

        return success, data

    if not check_grid_ids_l0_index(grid_ids):
        success = False
        data = ['Grid IDs contain invalid indices']

        return success, data

    return True, None

'''
Funtions for geographic transformations on to, and off of the EASE-DGGS.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''
from geopandas import GeoSeries
from shapely.geometry import Point

from vgrid.utils.easedggs.constants import grid_spec, levels_specs, ease_crs, geo_crs, cell_scale_factors

from vgrid.utils.easedggs.dggs.utils import shift_range_ease, shift_range_grid_xy

def ease_coord_to_grid_xy_coord(x_ease, y_ease):
    '''
    Convert DGGS cell coordinates to EASE Grid v2 coordinates

    Parameters
    ----------
    x : float
       EASE grid v2 x coordinate.
    y : float
        EASE grid v2 y coordinate.

    Returns
    -------
    (x_grid, y_grid)
    x_grid : float
        The transformed EASE Grid x cooridinate, in grid space units
    y_grid :  float
        The transformed EASE Grid x cooridinate, in grid space units
    '''
    if not isinstance(x_ease, float) and not isinstance(y_ease, float):
        return False

    return (shift_range_ease(x_ease, 'x'), shift_range_ease(y_ease, 'y'))

def grid_xy_coord_to_ease_coord(x_grid, y_grid):
    '''
    Convert EASE Grid v2 coordinates to DGGS cell

    Parameters
    ----------
    x_grid : float
       EASE-DGGS x coordinate.
    y_grid : float
        EASE-DGGS y coordinate.

    Returns
    -------
    (x_ease, y_ease) : tuple
        The transformed EASE-DGGS x,y cooridinate in as EASE v2 x,y coordinates

    '''
    if not isinstance(x_grid, float) and not isinstance(y_grid, float):
        return False

    return (shift_range_grid_xy(x_grid, 'x'), shift_range_grid_xy(y_grid, 'y'))


def coords_ease_to_coords_grid(coords_ease, level = 0):
    '''
    Convert GeoSeries of EASE Grid v2 coordinates to EASE-DGGS coordinates

    Parameters
    ----------
    coords_ease : Geopandas GeoSeries
       GeoSeries of EASE Grid coordinate pair points (x y).
    level : str
        EASE-DGGS level.

    Returns
    -------
    coords_grid: Geopandas GeoSeries
        GeoSeries of transformed coodinante pair points in the EASE-DGGS coordinate system.
    '''
    if not isinstance(coords_ease, GeoSeries):
        return False

    coords_grid = coords_ease.apply(lambda coord: Point(shift_range_ease(coord.x, 'x'), shift_range_ease(coord.y, 'y')))

    return coords_grid

def coords_lon_lat_to_coords_ease(coords_lon_lat, source_crs = geo_crs, target_crs = ease_crs):
    '''
    Convience function for converting between geographic and EASE grid cooridnates.
    Parameters
    ----------
    coords_lon_lat : list
       List of coordinate pairs (longitude, latitude) of geographic cooridnates.
    source_crs : int
        The EPSG code for the source geographic coordinate pair. Default is 4326.
    target_crs : int
        The EPSG code for the target (EASE v2) coordinate pair. Default is ease_crs.

    Returns
    -------
    coords_ease : list
        List of transformed coordinate pairs (lon, lat), now in EASE Grid v2 cooridnates(ease_x, ease_y).
    '''

    if not isinstance(coords_lon_lat, list):
        return False

    coords_lon_lat = GeoSeries([Point(coord[0], coord[1]) for coord in coords_lon_lat], crs = source_crs)
    coords_ease = coords_lon_lat.to_crs(target_crs)

    return coords_ease

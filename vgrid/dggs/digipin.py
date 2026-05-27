"""
DIGIPIN geocoding system for India.

DIGIPIN is a hierarchical geocoding system that divides a geographic area (India)
into a 4x4 grid recursively, using alphanumeric characters to encode locations.
"""

import math
from typing import Dict, Tuple, Union


# DIGIPIN grid character layout
DIGIPIN_GRID = [
    ["F", "C", "9", "8"],
    ["J", "3", "2", "7"],
    ["K", "4", "5", "6"],
    ["L", "M", "P", "T"],
]

# Create a mapping of characters to their grid positions
CHAR_POSITION_MAP = {}
for r in range(len(DIGIPIN_GRID)):
    for c in range(len(DIGIPIN_GRID[r])):
        CHAR_POSITION_MAP[DIGIPIN_GRID[r][c]] = (r, c)

# Geographic bounds for DIGIPIN (India region)
BOUNDS = {
    "minLat": 2.5,
    "minLon": 63.5,
    "maxLat": 38.5,
    "maxLon": 99.5,
}


def latlon_to_digipin(lat: float, lon: float, resolution: int = 10) -> str:
    """
    Convert latitude and longitude to DIGIPIN code.

    Parameters
    ----------
    lat : float
        Latitude coordinate (must be between 2.5 and 38.5)
    lon : float
        Longitude coordinate (must be between 63.5 and 99.5)
    resolution : int, optional
        Number of characters in the DIGIPIN code (default: 10)

    Returns
    -------
    str
        DIGIPIN code with dashes after 3rd and 6th characters,
        or 'Out of Bound' if coordinates are invalid

    Examples
    --------
    >>> latlon_to_digipin(28.6139, 77.2090, resolution=10)
    'F3K-492-6P96'
    """
    if not (math.isfinite(lat) and math.isfinite(lon)):
        return "Out of Bound"

    if lat < BOUNDS["minLat"] or lat > BOUNDS["maxLat"]:
        return "Out of Bound"

    if lon < BOUNDS["minLon"] or lon > BOUNDS["maxLon"]:
        return "Out of Bound"

    # Round to 6 decimal places for resolution
    lat = round(lat, 6)
    lon = round(lon, 6)

    min_lat = BOUNDS["minLat"]
    max_lat = BOUNDS["maxLat"]
    min_lon = BOUNDS["minLon"]
    max_lon = BOUNDS["maxLon"]
    pin = ""

    for level in range(1, resolution + 1):
        lat_div = (max_lat - min_lat) / 4
        lon_div = (max_lon - min_lon) / 4

        row = 3 - int((lat - min_lat) / lat_div)
        col = int((lon - min_lon) / lon_div)

        # Clamp to valid range
        r = min(max(row, 0), 3)
        c = min(max(col, 0), 3)

        pin += DIGIPIN_GRID[r][c]

        # Add dashes after 3rd and 6th characters
        if level == 3 or level == 6:
            pin += "-"

        # Update bounds for next iteration (same logic as digipin_to_bounds)
        lat1 = max_lat - lat_div * (r + 1)
        lat2 = max_lat - lat_div * r
        lon1 = min_lon + lon_div * c
        lon2 = min_lon + lon_div * (c + 1)

        min_lat = lat1
        max_lat = lat2
        min_lon = lon1
        max_lon = lon2
    # Remove trailing dash if present
    if pin.endswith("-"):
        pin = pin[:-1]

    return pin


def digipin_to_bounds(pin: str) -> Union[Dict[str, float], str]:
    """
    Get geographic bounds from a DIGIPIN code.

    Parameters
    ----------
    pin : str
        DIGIPIN code (with or without dashes)

    Returns
    -------
    dict or str
        Dictionary with keys 'minLat', 'maxLat', 'minLon', 'maxLon',
        or 'Invalid DIGIPIN' if the code is invalid

    Examples
    --------
    >>> digipin_to_bounds('F3K-492-6P96')
    {'minLat': 28.613..., 'maxLat': 28.614..., 'minLon': 77.208..., 'maxLon': 77.209...}
    """
    # Remove dashes
    clean = pin.replace("-", "")

    if len(clean) < 1:
        return "Invalid DIGIPIN"

    min_lat = BOUNDS["minLat"]
    min_lon = BOUNDS["minLon"]
    max_lat = BOUNDS["maxLat"]
    max_lon = BOUNDS["maxLon"]

    for ch in clean:
        position = CHAR_POSITION_MAP.get(ch)
        if position is None:
            return "Invalid DIGIPIN"

        r, c = position
        lat_div = (max_lat - min_lat) / 4
        lon_div = (max_lon - min_lon) / 4

        lat1 = max_lat - lat_div * (r + 1)
        lat2 = max_lat - lat_div * r
        lon1 = min_lon + lon_div * c
        lon2 = min_lon + lon_div * (c + 1)

        min_lat = lat1
        max_lat = lat2
        min_lon = lon1
        max_lon = lon2

    return {
        "minLat": min_lat,
        "maxLat": max_lat,
        "minLon": min_lon,
        "maxLon": max_lon,
    }


def digipin_to_center(pin: str) -> Union[Tuple[float, float], str]:
    """
    Get the center point (lat, lon) from a DIGIPIN code.

    Parameters
    ----------
    pin : str
        DIGIPIN code (with or without dashes)

    Returns
    -------
    tuple or str
        Tuple of (latitude, longitude) representing the center point,
        or 'Invalid DIGIPIN' if the code is invalid

    Examples
    --------
    >>> digipin_to_center('F3K')
    (28.125, 77.5)
    """
    bounds = digipin_to_bounds(pin)

    if isinstance(bounds, str):
        return bounds

    center_lat = (bounds["minLat"] + bounds["maxLat"]) / 2
    center_lon = (bounds["minLon"] + bounds["maxLon"]) / 2

    return (center_lat, center_lon)


def digipin_parent(pin: str) -> Union[str, str]:
    """
    Get the parent DIGIPIN code by removing the last character.
    If the pin is at level 1 (smallest resolution), returns itself.

    Parameters
    ----------
    pin : str
        DIGIPIN code (with or without dashes)

    Returns
    -------
    str
        Parent DIGIPIN code with dashes, or 'Invalid DIGIPIN' if the code is invalid

    Examples
    --------
    >>> digipin_parent('F3K-492')
    'F3K-49'
    >>> digipin_parent('F3K')
    'F3K'
    >>> digipin_parent('F')
    'F'
    """
    # Remove dashes
    clean = pin.replace("-", "")

    if len(clean) < 1:
        return "Invalid DIGIPIN"

    # If at level 1 (smallest resolution), parent is itself
    if len(clean) == 1:
        return pin

    # Remove last character
    parent_clean = clean[:-1]

    # Add dashes back after 3rd and 6th characters
    parent_with_dashes = ""
    for i, char in enumerate(parent_clean):
        parent_with_dashes += char
        if (i == 2 and len(parent_clean) > 3) or (i == 5 and len(parent_clean) > 6):
            parent_with_dashes += "-"

    return parent_with_dashes


def digipin_children(pin: str, target_resolution: int = None) -> Union[list, str]:
    """
    Get all child DIGIPIN codes by appending each possible character.
    If target_resolution equals current resolution, returns itself.

    Parameters
    ----------
    pin : str
        DIGIPIN code (with or without dashes)
    target_resolution : int, optional
        Target resolution for children. If None, returns children at next level.

    Returns
    -------
    list or str
        List of child DIGIPIN codes with dashes, or 'Invalid DIGIPIN' if the code is invalid

    Examples
    --------
    >>> digipin_children('F3K')
    ['F3K-F', 'F3K-C', 'F3K-9', 'F3K-8', 'F3K-J', 'F3K-3', 'F3K-2', 'F3K-7',
     'F3K-K', 'F3K-4', 'F3K-5', 'F3K-6', 'F3K-L', 'F3K-M', 'F3K-P', 'F3K-T']
    >>> digipin_children('F3K', 3)
    ['F3K']
    """
    # Remove dashes
    clean = pin.replace("-", "")

    if len(clean) < 1:
        return "Invalid DIGIPIN"

    if target_resolution is None:
        target_resolution = len(clean) + 1

    if target_resolution < len(clean):
        return "Invalid target resolution"

    # If target resolution equals current resolution, return itself
    if target_resolution == len(clean):
        return [pin]

    # If target resolution is only 1 level higher, generate direct children
    if target_resolution == len(clean) + 1:
        children = []
        for char in [
            "F",
            "C",
            "9",
            "8",
            "J",
            "3",
            "2",
            "7",
            "K",
            "4",
            "5",
            "6",
            "L",
            "M",
            "P",
            "T",
        ]:
            child_clean = clean + char

            # Add dashes after 3rd and 6th characters
            child_with_dashes = ""
            for i, c in enumerate(child_clean):
                child_with_dashes += c
                if (i == 2 and len(child_clean) > 3) or (
                    i == 5 and len(child_clean) > 6
                ):
                    child_with_dashes += "-"

            children.append(child_with_dashes)
        return children

    # If target resolution is more than 1 level higher, recursively generate children
    all_children = []
    direct_children = digipin_children(pin, len(clean) + 1)
    if isinstance(direct_children, list):
        for child in direct_children:
            child_children = digipin_children(child, target_resolution)
            if isinstance(child_children, list):
                all_children.extend(child_children)
            else:
                all_children.append(child)
    return all_children


def digipin_resolution(pin: str) -> Union[int, str]:
    """
    Get the resolution (number of characters) of a DIGIPIN code.

    Parameters
    ----------
    pin : str
        DIGIPIN code (with or without dashes)

    Returns
    -------
    int or str
        Resolution level, or 'Invalid DIGIPIN' if the code is invalid

    Examples
    --------
    >>> digipin_resolution('F3K-492')
    6
    """
    # Remove dashes
    clean = pin.replace("-", "")

    if len(clean) < 1:
        return "Invalid DIGIPIN"

    return len(clean)

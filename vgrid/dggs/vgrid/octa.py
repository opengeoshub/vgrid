import json
import math
from dataclasses import dataclass
from typing import List, Optional
from vgrid.utils.constants import AUTHALIC_RADIUS
# Reference: https://www.benthamopenarchives.com/contents/pdf/TOCSJ/TOCSJ-8-889.pdf
# An Equal Area Quad-tree Global Discrete Grid Subdivision Algorithm Based on Latitude and Longitude Lines
@dataclass
class Cell:
    lambda1: float  # λ₁: longitude boundary 1 (from parent)
    lambda2: float  # λ₂: longitude boundary 2 (from parent)
    phi1: float     # φ₁: base latitude or pole latitude (from parent)
    phi2: float     # φ₂: first subdivision latitude (from parent)
    phi3: float     # φ₃: second subdivision latitude (from parent)
    quad_code: str
    is_triangle: bool = False  # True if this is a triangular cell at the pole
    pole_lat: Optional[float] = None  # Latitude of the pole (90 or -90) if triangular


# Earth radius in meters (WGS84)
# EARTH_RADIUS = 6378137.0
EARTH_RADIUS = AUTHALIC_RADIUS
def calculate_cell_area_at_level(level: int, max_level: int) -> float:
    """
    Calculate the area F of a cell at a given subdivision level.
    
    Since we're doing equal-area subdivision, the area at level L is:
    Total sphere area / (8 * 4^L)
    
    Total sphere area = 4 * π * R²
    """
    total_sphere_area = 4 * math.pi * EARTH_RADIUS * EARTH_RADIUS
    num_cells_at_level = 8 * (4 ** level)
    return total_sphere_area / num_cells_at_level


def calculate_trapezoid_subdivision_latitude(
    lat1: float,
    lon1: float,
    lon2: float,
    level: int
) -> float:
    """
    Calculate subdivision latitude line for a spherical trapezoid cell.
    
    According to the paper (Section 3.2, part 2):
    - λ3 = (λ1 + λ2) / 2 (mean longitude) - handled separately
    - φ3 = arcsin(sin(φ1) + F/(R²(λ2-λ1)))
    
    Where:
    - F is the area of the divided trapezoid cell (child cell at next level)
    - R is the radius of the sphere
    - λ1, λ2 are the longitude lines of the trapezoid
    - φ1 is the known latitude line (lat_min for the cell)
    
    Args:
        lat1: Known latitude φ1 (typically lat_min)
        lon1: First longitude boundary λ1
        lon2: Second longitude boundary λ2
        level: Current subdivision level
        
    Returns:
        lat3 - the subdivision latitude line φ3
    """
    # Calculate area F for a child cell at the next level
    F = calculate_cell_area_at_level(level + 1, level + 1)
    
    # Convert longitudes to radians and get difference
    lon_diff_rad = abs(math.radians(lon2 - lon1))
    if lon_diff_rad > math.pi:
        lon_diff_rad = 2 * math.pi - lon_diff_rad
    
    # Avoid division by zero
    if lon_diff_rad < 1e-10:
        lon_diff_rad = math.pi / 2.0  # Default to 90 degrees
    
    # Convert latitude to radians
    lat1_rad = math.radians(lat1)
    sin_phi1 = math.sin(lat1_rad)
    
    # Calculate φ3: arcsin(sin(φ1) + F/(R²(λ2-λ1)))
    denominator = EARTH_RADIUS * EARTH_RADIUS * lon_diff_rad
    if denominator < 1e-10:
        denominator = 1e-10
    
    term = F / denominator
    sin_phi3 = sin_phi1 + term
    
    # Clamp to valid range [-1, 1]
    sin_phi3 = max(-1.0, min(1.0, sin_phi3))
    lat3 = math.degrees(math.asin(sin_phi3))
    
    return lat3


def calculate_triangle_subdivision_latitudes(
    lat_pole: float, 
    lon1: float, 
    lon2: float, 
    level: int
) -> tuple[float, float]:
    """
    Calculate subdivision latitude lines for a spherical triangle cell.
    
    According to the paper (Section 3.2):
    - λ3 = (λ1 + λ2) / 2 (mean longitude)
    - Area formula: F/2 = R²(λ₂ - λ₁)(sinφ₂ - sinφ₁)
    - Rearranging: φ2 = arcsin(sin(φ1) + F/(2R²(λ2-λ1)))
    - φ3 = arcsin [ F / (2R²(λ₂ - λ₁)) + sin(φ₁) ]
    
    Where:
    - F is the area of the divided triangle cell (child cell at next level)
    - R is the radius of the sphere
    - λ1, λ2 are the longitude lines of the triangle
    - φ1 is the latitude line of the triangle (at the pole, 90 or -90)
    
    Args:
        lat_pole: Latitude of the pole (90 or -90)
        lon1: First longitude boundary
        lon2: Second longitude boundary
        level: Current subdivision level
        
    Returns:
        Tuple of (lat2, lat3) - the two subdivision latitude lines
        For north pole: lat2 and lat3 are between pole (90) and base
        For south pole: lat2 and lat3 are between base and pole (-90)
    """
    # Calculate area F for a child cell at the next level
    F = calculate_cell_area_at_level(level + 1, level + 1)
    
    # Convert longitudes to radians and get difference
    lon_diff_rad = abs(math.radians(lon2 - lon1))
    if lon_diff_rad > math.pi:
        lon_diff_rad = 2 * math.pi - lon_diff_rad
    
    # Avoid division by zero
    if lon_diff_rad < 1e-10:
        lon_diff_rad = math.pi / 2.0  # Default to 90 degrees
    
    # Convert pole latitude to radians
    lat_pole_rad = math.radians(lat_pole)
    sin_phi1 = math.sin(lat_pole_rad)  # sin(90) = 1, sin(-90) = -1
    
    # Calculate φ2 using the area formula: F/2 = R²(λ₂ - λ₁)(sinφ₂ - sinφ₁)
    # Rearranging: sinφ₂ = sinφ₁ + F / (2R²(λ₂ - λ₁))
    # So: φ₂ = arcsin(sinφ₁ + F/(2R²(λ₂-λ₁)))
    denominator = EARTH_RADIUS * EARTH_RADIUS * lon_diff_rad
    if denominator < 1e-10:
        denominator = 1e-10
    
    term1 = F / (2 * denominator)  # F / (2 * R² * (λ₂ - λ₁))
    sin_phi2 = sin_phi1 + term1
    
    # Clamp to valid range [-1, 1]
    sin_phi2 = max(-1.0, min(1.0, sin_phi2))
    lat2 = math.degrees(math.asin(sin_phi2))
    
    # Calculate φ3: arcsin [ F / (2R²(λ₂ - λ₁)) + sin(φ₁) ]
    # Formula: φ₃ = arcsin [ F / (2R²(λ₂ - λ₁)) + sin(φ₁) ]
    term2 = F / (2 * denominator)  # F / (2 * R² * (λ₂ - λ₁))
    sin_phi3 = sin_phi1 + term2
    
    # Clamp to valid range
    sin_phi3 = max(-1.0, min(1.0, sin_phi3))
    lat3 = math.degrees(math.asin(sin_phi3))
    
    return lat2, lat3


def get_cell_boundaries(cell: Cell) -> tuple[float, float, float, float]:
    """
    Calculate cell boundaries (lat_min, lat_max, lon_min, lon_max) from lambda1, lambda2, phi1, phi2, phi3.
    
    Based on the child code (last character of quad_code) and whether parent was triangle pattern:
    - Triangle pattern (from octant or triangle cell):
      - "0": Triangle at pole (phi3 to pole)
      - "1": Middle trapezoid (phi2 to phi3)
      - "2": Left bottom/top trapezoid (phi1 to phi2, lambda1 to lambda3)
      - "3": Right bottom/top trapezoid (phi1 to phi2, lambda3 to lambda2)
    - Trapezoid subdivision (4 children):
      - "0": SW (phi1 to phi3, lambda1 to lambda3)
      - "1": SE (phi1 to phi3, lambda3 to lambda2)
      - "2": NW (phi3 to phi2, lambda1 to lambda3)
      - "3": NE (phi3 to phi2, lambda3 to lambda2)
    """
    lambda3 = (cell.lambda1 + cell.lambda2) / 2.0  # Subdivision longitude
    
    if len(cell.quad_code) <= 1:
        # Base octant (level 0): use phi1 as base (0 for octants 0-7)
        octant = int(cell.quad_code[0]) if cell.quad_code else 0
        
        if octant < 4:
            # Northern hemisphere octant: phi1 (0) to pole (90)
            return (cell.phi1, 90.0, cell.lambda1, cell.lambda2)
        else:
            # Southern hemisphere octant: pole (-90) to phi1 (0)
            return (-90.0, cell.phi1, cell.lambda1, cell.lambda2)
    
    child_code = cell.quad_code[-1]
    
    # Check if this cell itself is a triangle
    if cell.is_triangle:
        # Triangle cell: from phi3 to pole
        if cell.pole_lat and cell.pole_lat > 0:
            # North pole
            return (cell.phi3, cell.pole_lat, cell.lambda1, cell.lambda2)
        else:
            # South pole
            pole = cell.pole_lat if cell.pole_lat else -90.0
            return (pole, cell.phi3, cell.lambda1, cell.lambda2)
    
    # Determine pattern by checking if parent was triangle pattern
    # Triangle pattern: parent has codes 0,1,2,3 where 0 is triangle, 1 is middle trapezoid
    # We can detect this by checking if parent cell was a triangle or if we're at level 1
    level = len(cell.quad_code) - 1
    is_triangle_pattern = False
    
    if level == 1:
        # Level 1 is always triangle pattern (from octant)
        is_triangle_pattern = True
    elif len(cell.quad_code) >= 2:
        # Check if parent was triangle (code ends in "0" at previous level)
        parent_code = cell.quad_code[:-1]
        if len(parent_code) > 0 and parent_code[-1] == "0":
            # Parent was triangle, so this is triangle pattern
            is_triangle_pattern = True
    
    if is_triangle_pattern:
        # Triangle pattern children
        # Determine hemisphere from octant number
        octant = int(cell.quad_code[0]) if cell.quad_code else 0
        is_north = octant < 4
        
        if child_code == "0":
            # Triangle at pole
            if is_north:
                pole = cell.pole_lat if cell.pole_lat else 90.0
                return (cell.phi3, pole, cell.lambda1, cell.lambda2)
            else:
                pole = cell.pole_lat if cell.pole_lat else -90.0
                return (pole, cell.phi3, cell.lambda1, cell.lambda2)
        elif child_code == "1":
            # Middle trapezoid
            if is_north:
                return (cell.phi2, cell.phi3, cell.lambda1, cell.lambda2)
            else:
                # For south: phi3 is closer to pole, phi2 is closer to equator
                return (cell.phi3, cell.phi2, cell.lambda1, cell.lambda2)
        elif child_code == "2":
            # Left bottom/top trapezoid
            if is_north:
                return (cell.phi1, cell.phi2, cell.lambda1, lambda3)
            else:
                # For south: phi1 is equator (0), phi2 is closer to pole
                return (cell.phi1, cell.phi2, cell.lambda1, lambda3)
        elif child_code == "3":
            # Right bottom/top trapezoid
            if is_north:
                return (cell.phi1, cell.phi2, lambda3, cell.lambda2)
            else:
                # For south: phi1 is equator (0), phi2 is closer to pole
                return (cell.phi1, cell.phi2, lambda3, cell.lambda2)
    
    # Trapezoid subdivision (4 equal parts)
    # phi3 is the subdivision latitude calculated from phi1
    if child_code == "0":
        # SW
        return (cell.phi1, cell.phi3, cell.lambda1, lambda3)
    elif child_code == "1":
        # SE
        return (cell.phi1, cell.phi3, lambda3, cell.lambda2)
    elif child_code == "2":
        # NW
        return (cell.phi3, cell.phi2, cell.lambda1, lambda3)
    elif child_code == "3":
        # NE
        return (cell.phi3, cell.phi2, lambda3, cell.lambda2)
    
    # Fallback
    return (cell.phi1, cell.phi2, cell.lambda1, cell.lambda2)


def is_triangular_cell(cell: Cell) -> bool:
    """
    Check if a cell is a triangular cell (touching a pole).
    
    According to the paper, triangular cells occur at the poles.
    """
    if cell.is_triangle:
        return True
    
    # Check boundaries to see if it touches a pole
    lat_min, lat_max, lon_min, lon_max = get_cell_boundaries(cell)
    touches_south_pole = abs(lat_min + 90.0) < 0.1 or lat_min <= -89.9
    touches_north_pole = abs(lat_max - 90.0) < 0.1 or lat_max >= 89.9
    
    if touches_south_pole and touches_north_pole:
        return False  # Cell spans both poles, not a triangle
    return touches_south_pole or touches_north_pole


def subdivide_triangle_pattern(cell: Cell, level: int, max_level: int) -> List[Cell]:
    """
    Subdivide a triangular cell using the pattern: 1 triangle + 1 trapezoid + 2 trapezoids.
    
    This pattern is used for:
    - First level octant subdivision
    - All subsequent triangular cell subdivisions
    
    Pattern:
    - 1 triangle at the pole (full width)
    - 1 trapezoid in the middle (full width)
    - 2 trapezoids at the base (split by longitude)
    
    Creates children with lambda1, lambda2, phi1, phi2, phi3 from parent.
    """
    # Check if we've exceeded max_level (subdivide up to and including max_level)
    if level > max_level:
        return [cell]
    
    # Get cell boundaries to determine pole and base
    lat_min, lat_max, lon_min, lon_max = get_cell_boundaries(cell)
    
    # Determine which pole we're at
    if cell.pole_lat:
        pole_lat = cell.pole_lat
        is_north_pole = pole_lat > 0
    else:
        # Check octant number to determine pole
        octant = int(cell.quad_code[0]) if cell.quad_code else 0
        if octant < 4:
            # Northern hemisphere octants (0-3)
            pole_lat = 90.0
            is_north_pole = True
        else:
            # Southern hemisphere octants (4-7)
            pole_lat = -90.0
            is_north_pole = False
    
    # Calculate longitude midpoint: λ3 = (λ1 + λ2) / 2
    lambda3 = (cell.lambda1 + cell.lambda2) / 2.0
    
    # Check if this is an existing triangle cell (not first-level octant)
    is_existing_triangle = cell.is_triangle
    
    # Get the base latitude (furthest from pole) - this is phi1
    phi1_base = cell.phi1
    
    # Calculate phi2 and phi3 for children based on parent's phi1
    # According to the user: for resolution 1, F should be area of the whole globe / 8 (octant area)
    # In general, F is the area of the parent cell being subdivided
    # When subdividing at level L, F is the area at level L-1 (the parent's level)
    if level == 1:
        # For resolution 1, F is the octant area (level 0)
        F = calculate_cell_area_at_level(0, 0)
    else:
        # For other levels, F is the parent cell's area (level L-1)
        F = calculate_cell_area_at_level(level - 1, level - 1)
    
    # Convert longitudes to radians and get difference
    lon_diff_rad = abs(math.radians(cell.lambda2 - cell.lambda1))
    if lon_diff_rad > math.pi:
        lon_diff_rad = 2 * math.pi - lon_diff_rad
    
    if lon_diff_rad < 1e-10:
        lon_diff_rad = math.pi / 2.0
    
    denominator = EARTH_RADIUS * EARTH_RADIUS * lon_diff_rad
    if denominator < 1e-10:
        denominator = 1e-10
    
    if is_existing_triangle:
        # For existing triangle cells, phi1 is the pole latitude
        # Calculate phi2 and phi3 from pole using the exact formulas from the paper:
        # φ₂ = arcsin [ F / (2R²(λ₂-λ₁)) + sinφ₁ ]
        # φ₃ = arcsin [ F / (4R²(λ₂-λ₁)) + sinφ₂ ]
        if is_north_pole:
            sin_phi1 = math.sin(math.radians(pole_lat))  # sin(90) = 1
            # phi2 = arcsin [ F / (2R²(λ₂-λ₁)) + sinφ₁ ]
            term_phi2 = F / (2 * denominator)
            sin_phi2 = term_phi2 + sin_phi1  # F/(2R²(λ₂-λ₁)) + sin(φ₁)
            sin_phi2 = max(-1.0, min(1.0, sin_phi2))
            phi2 = math.degrees(math.asin(sin_phi2))
            
            # phi3 = arcsin [ F / (4R²(λ₂-λ₁)) + sinφ₂ ]
            term_phi3 = F / (4 * denominator)
            sin_phi3 = term_phi3 + sin_phi2  # F/(4R²(λ₂-λ₁)) + sin(φ₂)
            sin_phi3 = max(-1.0, min(1.0, sin_phi3))
            phi3 = math.degrees(math.asin(sin_phi3))
        else:
            sin_phi1 = math.sin(math.radians(pole_lat))  # sin(-90) = -1
            # phi2 = arcsin [ F / (2R²(λ₂-λ₁)) + sinφ₁ ]
            # For southern hemisphere, sin(φ₁) is negative, so we add the term
            term_phi2 = F / (2 * denominator)
            sin_phi2 = term_phi2 + sin_phi1  # F/(2R²(λ₂-λ₁)) + sin(φ₁)
            sin_phi2 = max(-1.0, min(1.0, sin_phi2))
            phi2 = math.degrees(math.asin(sin_phi2))
            
            # phi3 = arcsin [ F / (4R²(λ₂-λ₁)) + sinφ₂ ]
            term_phi3 = F / (4 * denominator)
            sin_phi3 = term_phi3 + sin_phi2  # F/(4R²(λ₂-λ₁)) + sin(φ₂)
            sin_phi3 = max(-1.0, min(1.0, sin_phi3))
            phi3 = math.degrees(math.asin(sin_phi3))
    else:
        # For first-level octant subdivision, phi1 is the base (equator = 0)
        # Using the exact formulas from the paper:
        # φ₂ = arcsin [ F / (2R²(λ₂-λ₁)) + sinφ₁ ]
        # φ₃ = arcsin [ F / (4R²(λ₂-λ₁)) + sinφ₂ ]
        if is_north_pole:
            sin_phi1 = math.sin(math.radians(phi1_base))  # sin(0) = 0
            # phi2 = arcsin [ F / (2R²(λ₂-λ₁)) + sinφ₁ ]
            term_phi2 = F / (2 * denominator)
            sin_phi2 = term_phi2 + sin_phi1  # F/(2R²(λ₂-λ₁)) + sin(φ₁)
            sin_phi2 = max(-1.0, min(1.0, sin_phi2))
            phi2 = math.degrees(math.asin(sin_phi2))
            
            # phi3 = arcsin [ F / (4R²(λ₂-λ₁)) + sinφ₂ ]
            term_phi3 = F / (4 * denominator)
            sin_phi3 = term_phi3 + sin_phi2  # F/(4R²(λ₂-λ₁)) + sin(φ₂)
            sin_phi3 = max(-1.0, min(1.0, sin_phi3))
            phi3 = math.degrees(math.asin(sin_phi3))
        else:
            # For south octants, phi1 is the base (equator = 0)
            # phi2 and phi3 should be negative (south of equator)
            sin_phi1 = math.sin(math.radians(phi1_base))  # sin(0) = 0
            # phi2 = arcsin [ F / (2R²(λ₂-λ₁)) + sinφ₁ ]
            # For south, we subtract instead of add to get negative values
            term_phi2 = F / (2 * denominator)
            sin_phi2 = sin_phi1 - term_phi2  # sin(φ₁) - F/(2R²(λ₂-λ₁)) for south
            sin_phi2 = max(-1.0, min(1.0, sin_phi2))
            phi2 = math.degrees(math.asin(sin_phi2))
            
            # phi3 = arcsin [ sin(φ₂) - F / (4R²(λ₂-λ₁)) ]
            term_phi3 = F / (4 * denominator)
            sin_phi3 = sin_phi2 - term_phi3  # sin(φ₂) - F/(4R²(λ₂-λ₁)) for south
            sin_phi3 = max(-1.0, min(1.0, sin_phi3))
            phi3 = math.degrees(math.asin(sin_phi3))
    
    # All children share the same lambda1, lambda2, phi1, phi2, phi3 from parent
    # Create children with these shared parameters
    if is_north_pole:
        children = [
            # Top triangle (full width, near pole) - code "0"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "0",
                 is_triangle=True, pole_lat=pole_lat),
            # Middle trapezoid (full width) - code "1"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "1"),
            # Left-bottom trapezoid (SW) - code "2"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "2"),
            # Right-bottom trapezoid (SE) - code "3"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "3"),
        ]
    else:
        children = [
            # Bottom triangle (full width, near pole) - code "0"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "0",
                 is_triangle=True, pole_lat=pole_lat),
            # Middle trapezoid (full width) - code "1"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "1"),
            # Left-top trapezoid (NW) - code "2"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "2"),
            # Right-top trapezoid (NE) - code "3"
            Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "3"),
        ]
    
    result = []
    for c in children:
        result.extend(subdivide_cell(c, level + 1, max_level))
    return result


def subdivide_triangular_cell(cell: Cell, level: int, max_level: int) -> List[Cell]:
    """
    Subdivide a triangular cell using the pattern: 1 triangle + 1 trapezoid + 2 trapezoids.
    
    This is the same pattern used at the first level, applied recursively to triangles.
    """
    return subdivide_triangle_pattern(cell, level, max_level)


def subdivide_octant_first_level(cell: Cell, max_level: int) -> List[Cell]:
    """
    Subdivide an octant at the first level (level 1) according to the paper.
    
    At resolution = 1, each octant is divided into:
    - 1 triangle at the top (near the pole)
    - 1 trapezoid in the middle
    - 2 trapezoids at the bottom (near equator)
    
    This creates 4 children total: 1 triangle + 3 trapezoids.
    """
    if max_level < 1:
        return [cell]
    
    # Use the same triangle pattern function, but start at level 1
    # The function will handle stopping at max_level
    return subdivide_triangle_pattern(cell, 1, max_level)


def subdivide_cell(cell: Cell, level: int, max_level: int) -> List[Cell]:
    """
    Recursively subdivide into 4 equal-area parts (quad-tree style).
    
    According to the paper, this maintains equal area cells at each subdivision level.
    
    Handles both triangular cells (at poles) and trapezoid (quadrilateral) cells.
    For trapezoid cells, uses the formula from Section 3.2:
    - Longitude division: λ3 = (λ1 + λ2) / 2 (mean longitude)
    - Latitude division: φ3 = arcsin(sin(φ1) + F/(R²(λ2-λ1)))
    """
    # Check if we've exceeded max_level (subdivide up to and including max_level)
    if level > max_level:
        return [cell]
    
    # Check if this is a triangular cell
    if is_triangular_cell(cell):
        return subdivide_triangular_cell(cell, level, max_level)
    
    # For trapezoid (quadrilateral) cells, use the paper's formula:
    # - Longitude division: λ3 = (λ1 + λ2) / 2 (mean longitude)
    # - Latitude division: φ3 = arcsin(sin(φ1) + F/(R²(λ2-λ1)))
    
    # Get cell boundaries to determine phi1 (base latitude)
    lat_min, lat_max, lon_min, lon_max = get_cell_boundaries(cell)
    
    # Calculate phi3 (subdivision latitude) from phi1
    # Use phi1 as the base latitude (lat_min for northern, lat_max for southern)
    phi1_base = cell.phi1
    
    # According to the paper, F is the area of the divided trapezoid cell (child cell at next level)
    # This ensures equal-area subdivision
    F = calculate_cell_area_at_level(level + 1, level + 1)
    
    # Convert longitudes to radians and get difference
    lon_diff_rad = abs(math.radians(cell.lambda2 - cell.lambda1))
    if lon_diff_rad > math.pi:
        lon_diff_rad = 2 * math.pi - lon_diff_rad
    
    if lon_diff_rad < 1e-10:
        lon_diff_rad = math.pi / 2.0
    
    denominator = EARTH_RADIUS * EARTH_RADIUS * lon_diff_rad
    if denominator < 1e-10:
        denominator = 1e-10
    
    # For trapezoid subdivision:
    # - phi1 is the base latitude (from parent)
    # - phi2 is the top boundary of the parent cell (should be cell.phi2 if available, otherwise calculate from boundaries)
    # - phi3 is the subdivision latitude calculated using the formula
    
    # phi2 should be the parent cell's top boundary
    # For a trapezoid cell, phi2 from the parent represents the top boundary
    if cell.phi2 != 0.0 or len(cell.quad_code) > 1:
        # Use parent's phi2 if it's set (non-zero or not base octant)
        phi2 = cell.phi2
    else:
        # Fallback: calculate from cell boundaries
        if lat_max > lat_min:
            phi2 = lat_max
        else:
            phi2 = lat_min
    
    # Calculate phi3: φ3 = arcsin(sin(φ1) + F/(R²(λ2-λ1)))
    # According to the paper, F is the area of the divided trapezoid cell (child at next level)
    sin_phi1 = math.sin(math.radians(phi1_base))
    term_phi3 = F / denominator
    sin_phi3 = sin_phi1 + term_phi3
    sin_phi3 = max(-1.0, min(1.0, sin_phi3))
    phi3 = math.degrees(math.asin(sin_phi3))
    
    # All children share the same lambda1, lambda2, phi1, phi2, phi3
    children = [
        Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "0"),  # SW
        Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "1"),  # SE
        Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "2"),  # NW
        Cell(cell.lambda1, cell.lambda2, cell.phi1, phi2, phi3, cell.quad_code + "3"),  # NE
    ]

    result = []
    for c in children:
        result.extend(subdivide_cell(c, level + 1, max_level))
    return result


def generate_equal_area_grid(max_level: int) -> List[Cell]:
    """
    Generate an equal-area quad-tree global discrete grid based on latitude and longitude lines.
    
    According to the paper, this algorithm:
    1. Starts with 8 base octants (hemispheres divided by equator and prime meridian)
    2. Recursively subdivides each cell into 4 equal-area children
    3. Maintains equal area at each subdivision level
    4. Has stable geometric distortion that converges to a ratio of ~2.0
    
    The subdivision uses equal-area latitude division rather than simple midpoint division
    to ensure all cells at the same level have equal area.
    """
    # Base octants: lambda1, lambda2, phi1 (0 for equator), phi2 (will be calculated), phi3 (will be calculated)
    # For base octants, phi2 and phi3 are not yet defined (set to 0 as placeholder)
    base_octants = [
        # Region 0: Lat 0° to +90°, Lon 0° to +90° (Northern Hemisphere, Eastern Hemisphere)
        Cell(0.0, 90.0, 0.0, 0.0, 0.0, "0"),
        # Region 1: Lat 0° to +90°, Lon +90° to +180°
        Cell(90.0, 180.0, 0.0, 0.0, 0.0, "1"),
        # Region 2: Lat 0° to +90°, Lon -180° to -90°
        Cell(-180.0, -90.0, 0.0, 0.0, 0.0, "2"),
        # Region 3: Lat 0° to +90°, Lon -90° to 0°
        Cell(-90.0, 0.0, 0.0, 0.0, 0.0, "3"),
        # Region 4: Lat -90° to 0°, Lon 0° to +90° (Southern Hemisphere, Eastern Hemisphere)
        Cell(0.0, 90.0, 0.0, 0.0, 0.0, "4"),
        # Region 5: Lat -90° to 0°, Lon +90° to +180°
        Cell(90.0, 180.0, 0.0, 0.0, 0.0, "5"),
        # Region 6: Lat -90° to 0°, Lon -180° to -90°
        Cell(-180.0, -90.0, 0.0, 0.0, 0.0, "6"),
        # Region 7: Lat -90° to 0°, Lon -90° to 0°
        Cell(-90.0, 0.0, 0.0, 0.0, 0.0, "7"),
    ]

    all_cells = []
    for cell in base_octants:
        # First level (level 1) has special subdivision pattern:
        # 1 triangle at pole + 1 trapezoid in middle + 2 trapezoids at equator
        all_cells.extend(subdivide_octant_first_level(cell, max_level))
    return all_cells


def calculate_cell_area(cell: Cell) -> float:
    """
    Calculate the area of a spherical cell.
    
    For a spherical rectangle: Area = R² * (λ₂ - λ₁) * (sin(φ₂) - sin(φ₁))
    """
    lat_min, lat_max, lon_min, lon_max = get_cell_boundaries(cell)
    
    lon_diff_rad = abs(math.radians(lon_max - lon_min))
    if lon_diff_rad > math.pi:
        lon_diff_rad = 2 * math.pi - lon_diff_rad
    
    lat_min_rad = math.radians(lat_min)
    lat_max_rad = math.radians(lat_max)
    
    area = EARTH_RADIUS * EARTH_RADIUS * lon_diff_rad * (math.sin(lat_max_rad) - math.sin(lat_min_rad))
    return area


def find_parent_cell(cell: Cell) -> Optional[Cell]:
    """
    Find the parent cell by regenerating the grid up to the parent level.
    
    Returns the parent Cell if found, None if cell is at level 0.
    """
    level = len(cell.quad_code) - 1 if len(cell.quad_code) > 0 else 0
    
    if level == 0:
        return None  # Base octant has no parent
    
    parent_level = level - 1
    parent_quad_code = cell.quad_code[:-1]  # Remove last character
    
    # Regenerate grid up to parent level to find the parent cell
    parent_cells = generate_equal_area_grid(parent_level)
    
    # Find the parent cell with matching quad_code
    for parent_cell in parent_cells:
        if parent_cell.quad_code == parent_quad_code:
            return parent_cell
    
    return None


def cell_to_geojson_feature(cell: Cell):
    """
    Convert a cell into a GeoJSON polygon feature.
    Polygon coordinates are in WGS84 (lat/lon).
    Includes F (area), lambda1, lambda2, phi1, phi2, phi3 in properties.
    """
    # Get cell boundaries from lambda1, lambda2, phi1, phi2, phi3
    lat_min, lat_max, lon_min, lon_max = get_cell_boundaries(cell)
    
    polygon = [
        [
            [lon_min, lat_min],
            [lon_max, lat_min],
            [lon_max, lat_max],
            [lon_min, lat_max],
            [lon_min, lat_min],
        ]
    ]
    
    # Calculate resolution from quad_code length
    resolution = len(cell.quad_code) - 1 if len(cell.quad_code) > 0 else 0
    
    # Calculate F: at resolution n, F is the area of its parent at resolution n-1
    # Since we're doing equal-area subdivision, F = area at resolution (n-1)
    if resolution == 0:
        # Base octant: F is the area of the entire octant (1/8 of sphere)
        F = calculate_cell_area(cell)
    else:
        # For cells at resolution n, F is the area of parent at resolution n-1
        F = calculate_cell_area_at_level(resolution - 1, resolution - 1)
    
    # Use cell's lambda1, lambda2, phi1, phi2, phi3 (from parent)
    properties = {
        "quad_code": cell.quad_code,
        "F": F,  # Area of the parent cell (m²)
        "lambda1": cell.lambda1,  # λ₁ (longitude boundary 1, degrees)
        "lambda2": cell.lambda2,  # λ₂ (longitude boundary 2, degrees)
        "phi1": cell.phi1,  # φ₁ (base latitude or pole latitude, degrees)
        "phi2": cell.phi2,  # φ₂ (subdivision latitude, degrees)
        "phi3": cell.phi3,  # φ₃ (subdivision latitude, degrees)
        "resolution": resolution,  # Subdivision resolution
    }

    return {
        "type": "Feature",
        "properties": properties,
        "geometry": {
            "type": "Polygon",
            "coordinates": polygon
        }
    }


def export_geojson(cells: List[Cell], filename="grid2_authalic.geojson"):
    features = [cell_to_geojson_feature(c) for c in cells]

    geojson = {
        "type": "FeatureCollection",
        "features": features
    }

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(geojson, f, indent=2)

    print(f"GeoJSON saved to {filename}")
    print(f"Total cells: {len(features)}")


if __name__ == "__main__":
    max_level = 2    # <-- change depth here
    cells = generate_equal_area_grid(max_level)
    export_geojson(cells)

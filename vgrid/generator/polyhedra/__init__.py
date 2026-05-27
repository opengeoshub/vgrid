"""
Polyhedra module for vgrid.

This module provides various polyhedra implementations used in discrete global grid systems (DGGS).
"""

from .cube import cube, cube_cli
from .octahedron import octahedron, octahedron_cli
from .tetrahedron import tetrahedron, tetrahedron_cli
from .dodecahedron import dodecahedron, dodecahedron_cli
from .fuller_icosahedron import fuller_icosahedron, fuller_icosahedron_cli
from .rhombic_icosahedron import rhombic_icosahedron, rhombic_icosahedron_cli

__all__ = [
    "cube",
    "cube_cli",
    "octahedron",
    "octahedron_cli",
    "tetrahedron",
    "dodecahedron",
    "dodecahedron_cli",
    "fuller_icosahedron",
    "fuller_icosahedron_cli",
    "rhombic_icosahedron",
    "rhombic_icosahedron_cli",
]

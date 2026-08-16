"""
DGGS Compact and Expand functions.

This submodule provides functions to compact and expand various discrete global grid systems (DGGS).
"""

from .h3compact import h3compact, h3expand, h3compact_cli, h3expand_cli
from .s2compact import s2compact, s2expand, s2compact_cli, s2expand_cli
from .a5compact import a5compact, a5expand, a5compact_cli, a5expand_cli
from .rhealpixcompact import (
    rhealpixcompact,
    rhealpixexpand,
    rhealpixcompact_cli,
    rhealpixexpand_cli,
)
from .isea4tcompact import (
    isea4tcompact,
    isea4texpand,
    isea4tcompact_cli,
    isea4texpand_cli,
)
from .isea3hcompact import (
    isea3hcompact,
    isea3hexpand,
    isea3hcompact_cli,
    isea3hexpand_cli,
)
from .easecompact import easecompact, easeexpand, easecompact_cli, easeexpand_cli
from .dggalcompact import dggalcompact, dggalexpand, dggalcompact_cli, dggalexpand_cli
from .qtmcompact import qtmcompact, qtmexpand, qtmcompact_cli, qtmexpand_cli
from .olccompact import olccompact, olcexpand, olccompact_cli, olcexpand_cli
from .geohashcompact import (
    geohashcompact,
    geohashexpand,
    geohashcompact_cli,
    geohashexpand_cli,
)
from .tilecodecompact import (
    tilecodecompact,
    tilecodeexpand,
    tilecodecompact_cli,
    tilecodeexpand_cli,
)
from .quadkeycompact import (
    quadkeycompact,
    quadkeyexpand,
    quadkeycompact_cli,
    quadkeyexpand_cli,
)
from .digipincompact import (
    digipincompact,
    digipinexpand,
    digipincompact_cli,
    digipinexpand_cli,
)
from .healpixcompact import (
    healpixcompact,
    healpixexpand,
    healpixcompact_cli,
    healpixexpand_cli,
)

__all__ = [
    "h3compact",
    "h3expand",
    "h3compact_cli",
    "h3expand_cli",
    "s2compact",
    "s2expand",
    "s2compact_cli",
    "s2expand_cli",
    "a5compact",
    "a5expand",
    "a5compact_cli",
    "a5expand_cli",
    "rhealpixcompact",
    "rhealpixexpand",
    "rhealpixcompact_cli",
    "rhealpixexpand_cli",
    "isea4tcompact",
    "isea4texpand",
    "isea4tcompact_cli",
    "isea4texpand_cli",
    "isea3hcompact",
    "isea3hexpand",
    "isea3hcompact_cli",
    "isea3hexpand_cli",
    "easecompact",
    "easeexpand",
    "easecompact_cli",
    "easeexpand_cli",
    "dggalcompact",
    "dggalexpand",
    "dggalcompact_cli",
    "dggalexpand_cli",
    "qtmcompact",
    "qtmexpand",
    "qtmcompact_cli",
    "qtmexpand_cli",
    "olccompact",
    "olcexpand",
    "olccompact_cli",
    "olcexpand_cli",
    "geohashcompact",
    "geohashexpand",
    "geohashcompact_cli",
    "geohashexpand_cli",
    "tilecodecompact",
    "tilecodeexpand",
    "tilecodecompact_cli",
    "tilecodeexpand_cli",
    "quadkeycompact",
    "quadkeyexpand",
    "quadkeycompact_cli",
    "quadkeyexpand_cli",
    "digipincompact",
    "digipinexpand",
    "digipincompact_cli",
    "digipinexpand_cli",
    "healpixcompact",
    "healpixexpand",
    "healpixcompact_cli",
    "healpixexpand_cli",
]

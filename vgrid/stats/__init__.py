"""
Statistics module for vgrid.

This module provides functions to calculate and display statistics for various
discrete global grid systems (DGGS), including cell counts, areas, and edge lengths.
"""

from .h3stats import h3stats, h3stats_cli, h3inspect, h3inspect_cli
from .s2stats import s2stats, s2stats_cli, s2inspect, s2inspect_cli
from .a5stats import a5stats, a5stats_cli, a5inspect, a5inspect_cli
from .rhealpixstats import (
    rhealpixstats,
    rhealpixstats_cli,
    rhealpixinspect,
    rhealpixinspect_cli,
)
from .isea4tstats import isea4tstats, isea4tstats_cli, isea4tinspect, isea4tinspect_cli
from .isea3hstats import isea3hstats, isea3hstats_cli, isea3hinspect, isea3hinspect_cli
from .easestats import easestats, easestats_cli, easeinspect, easeinspect_cli
from .qtmstats import qtmstats, qtmstats_cli, qtminspect, qtminspect_cli
from .olcstats import olcstats, olcstats_cli, olcinspect, olcinspect_cli
from .geohashstats import (
    geohashstats,
    geohashstats_cli,
    geohashinspect,
    geohashinspect_cli,
)
from .georefstats import georefstats, georefstats_cli, georefinspect, georefinspect_cli
from .mgrsstats import mgrsstats, mgrsstats_cli
from .tilecodestats import (
    tilecodestats,
    tilecodestats_cli,
    tilecodeinspect,
    tilecodeinspect_cli,
)
from .quadkeystats import (
    quadkeystats,
    quadkeystats_cli,
    quadkeyinspect,
    quadkeyinspect_cli,
)
from .maidenheadstats import (
    maidenheadstats,
    maidenheadstats_cli,
    maidenheadinspect,
    maidenheadinspect_cli,
)
from .garsstats import garsstats, garsstats_cli, garsinspect, garsinspect_cli
from .digipinstats import (
    digipinstats,
    digipinstats_cli,
    digipininspect,
    digipininspect_cli,
)
from .dggalstats import dggalstats, dggalstats_cli, dggalinspect, dggalinspect_cli
from .dggridstats import dggridstats, dggridstats_cli, dggridinspect, dggridinspect_cli

__all__ = [
    "h3stats",
    "h3stats_cli",
    "h3inspect",
    "h3inspect_cli",
    "s2stats",
    "s2stats_cli",
    "s2inspect",
    "s2inspect_cli",
    "a5stats",
    "a5stats_cli",
    "a5inspect",
    "a5inspect_cli",
    "rhealpixstats",
    "rhealpixstats_cli",
    "rhealpixinspect",
    "rhealpixinspect_cli",
    "isea4tstats",
    "isea4tstats_cli",
    "isea4tinspect",
    "isea4tinspect_cli",
    "isea3hstats",
    "isea3hstats_cli",
    "isea3hinspect",
    "isea3hinspect_cli",
    "easestats",
    "easestats_cli",
    "easeinspect",
    "easeinspect_cli",
    "qtmstats",
    "qtmstats_cli",
    "qtminspect",
    "qtminspect_cli",
    "olcstats",
    "olcstats_cli",
    "olcinspect",
    "olcinspect_cli",
    "geohashstats",
    "geohashstats_cli",
    "geohashinspect",
    "geohashinspect_cli",
    "georefstats",
    "georefstats_cli",
    "georefinspect",
    "georefinspect_cli",
    "mgrsstats",
    "mgrsstats_cli",
    "tilecodestats",
    "tilecodestats_cli",
    "tilecodeinspect",
    "tilecodeinspect_cli",
    "quadkeystats",
    "quadkeystats_cli",
    "quadkeyinspect",
    "quadkeyinspect_cli",
    "maidenheadstats",
    "maidenheadstats_cli",
    "maidenheadinspect",
    "maidenheadinspect_cli",
    "garsstats",
    "garsstats_cli",
    "garsinspect",
    "garsinspect_cli",
    "digipinstats",
    "digipinstats_cli",
    "digipininspect",
    "digipininspect_cli",
    "dggalstats",
    "dggalstats_cli",
    "dggalinspect",
    "dggalinspect_cli",
    "dggridstats",
    "dggridstats_cli",
    "dggridinspect",
    "dggridinspect_cli",
]

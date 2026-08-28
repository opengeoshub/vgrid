"""pytest session header: tolerance table."""

from _common import format_tolerance_table


def pytest_report_header(config):
    return ["validation tolerances", format_tolerance_table()]

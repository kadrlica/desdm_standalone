#!/usr/bin/env python
"""
DESDM standalone package
"""
__author__ = "Alex Drlica-Wagner"
__email__ = "kadrlica@fnal.gov"

try:
    from .version import __version__
except ImportError:
    from .get_version import get_version
    __version__ = get_version()


#!/usr/bin/env python
"""
DESDM standalone package
"""
__author__ = "Alex Drlica-Wagner"
__email__ = "kadrlica@fnal.gov"

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

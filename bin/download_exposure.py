#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

from standalone import noao

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('expnum',type=int)
    args = parser.parse_args()

    noao.download_exposure(expnum)

#!/usr/bin/env python
"""
Module for dealing with observing epochs.

Some documentation:
Official DES Observing Epochs:
https://opensource.ncsa.illinois.edu/confluence/display/DESDM/Y2A1+Observing+Epochs

DESDM Observing Events:
https://opensource.ncsa.illinois.edu/confluence/display/DESDM/Observing+Events

Operations Configuration Changes:
https://cdcvs.fnal.gov/redmine/projects/desops/wiki/Log_of_Configuration_Changes

The Dark Energy Survey and Operations: Year 1
http://adsabs.harvard.edu/abs/2014SPIE.9149E..0VD

"""
__author__ = "Alex Drlica-Wagner"

from collections import OrderedDict as odict
from dateutil.parser import parse as  dtparse

# Official dates can (mostly) be found in the documentation. Here we
# extend the epochs to be continuous across time and to extend into Y3.
#https://opensource.ncsa.illinois.edu/confluence/display/DESDM/Y2A1+Observing+Epochs

EPOCHS = odict([
        (1, dict(year='SV',epoch='E1',
                 tstart=dtparse('20120912'),tstop=dtparse('20121229'))),
        (2, dict(year='SV',epoch='E2',
                 tstart=dtparse('20121229'),tstop=dtparse('20130830'))),
        (3, dict(year='Y1',epoch='E1',
                 tstart=dtparse('20130831'),tstop=dtparse('20131129'))),
        (4, dict(year='Y1',epoch='E2',
                 tstart=dtparse('20131129'),tstop=dtparse('20140815'))),
        (5, dict(year='Y2',epoch='E1',
                 tstart=dtparse('20140815'),tstop=dtparse('20141130'))),
        # Currently Y2_E2 extends through Y3...
        (6, dict(year='Y2',epoch='E2',
                 tstart=dtparse('20141130'),tstop=dtparse('20160801'))),
])

def nite2key(nite):
    if isinstance(nite,basestring):
        nite = dtparse(nite)

    for k,v in EPOCHS.items():
        if v['tstart'] < nite < v['tstop']:
            return k
    msg = "Nite outside of EPOCHS: %s"%nite
    raise Exception(msg)

def nite2val(nite):
    return EPOCHS[nite2iter(nite)]

def nite2year(nite):
    return nite2val(nite)['year']

def nite2epoch(nite):
    return nite2val['epoch']


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

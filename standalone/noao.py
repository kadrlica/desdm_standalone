#!/usr/bin/env python
"""
Interface to NOAO data archive.

Query URL: http://portal-nvo.noao.edu/search/query
Schema Doc: http://portal-nvo.noao.edu/search/show_schema

Step by step instructions: 
1) Create the query with 'get_noao_query'
2) Navigate to 'http://portal-nvo.noao.edu/search/query', and enter
the query in the 'Advanced Query Form' tab.
3) Click 'Download rows as VOTable' and select 'All rows'.
"""
__author__ = "Alex Drlica-Wagner"

from datetime import date
from templates import noao_query,noao_curl
import subprocess

import numpy as np
import numpy.lib.recfunctions as recfuncs
import astropy.io.votable as vot

NOAO_URL = 'http://portal-nvo.noao.edu/search/query'
VOTABLE = '/home/s1/kadrlica/projects/standalone/data/votables/votable_20160131.vot'

def get_noao_query(**kwargs):
    kwargs = get_noao_query_kwargs(**kwargs)
    return noao_query.format(**kwargs)

def get_noao_query_kwargs(**kwargs):
    """
    Get the NOAO download query.
    """
    # Some columns are required
    required = [
        'reference', 
        'release_date', 
        'start_date', 
        'filesize', 
        'dtpropid', 
        'md5sum'
    ]

    defaults = dict(today=date.today(),exptime=60,filters=('g','r','i','z','Y'),
                    limit=250000)
    defaults['columns'] = [
        'reference', 
        'dtpropid', 
        'release_date', 
        'start_date', 
        'date_obs', 
        'instrument', 
        'ra', 
        'dec', 
        'filter', 
        'exposure', 
        'obstype', 
        'proctype', 
        'dtacqnam AS original_file', 
        'reference AS archive_file',
        'filesize',
        ]

    for k,v in defaults.items():
        kwargs.setdefault(k,v)

    kwargs['columns'] = map(str.lower,kwargs['columns'])
    kwargs['columns'] += [c for c in required if c not in kwargs['columns']]

    if not isinstance(kwargs['columns'],basestring):
        kwargs['columns'] = ','.join(kwargs['columns'])
    if not isinstance(kwargs['filters'],basestring):
        kwargs['filters'] = ','.join(["'%s'"%f for f in kwargs['filters']])

    return kwargs

def get_noao_exposures(votable=None):
    pass

def create_expnum(data):
    col = 'original_file'
    dtype ='S%i'%(len(max(data[col], key=len)))
    filenames = data[col].data.astype(dtype)
    basenames = np.char.rpartition(filenames,'/')[:,-1]
    splitexts = np.char.strip(basenames,'.fits.fz')
    expnum = np.char.rpartition(splitexts,'_')[:,-1].astype(int)
    return expnum

def create_nite(data):
    """ This is a nitemare. """
    col = 'date_obs'
    datetime = np.array(data[col],dtype='datetime64')
    date,sep,time = np.char.partition(data[col].astype('S30'),' ').T
    nite = np.char.replace(date,'-','').astype(int)
    date = np.array(np.char.add(date,' 00:00:00.0'),dtype='datetime64')
    last = datetime - date < np.timedelta64(14,'h')
    nite[last] -= 1
    return nite

def vot2npy(votable):
    if isinstance(votable,basestring):
        data = vot.parse_single_table(votable).array
    else:
        data = votable

    expnum = create_expnum(data)
    nite = create_nite(data)
    out = np.empty(len(data),dtype=data.dtype.descr+[('expnum',int),('nite',int)])
    out[:] = data[:]
    out['expnum'] = expnum
    out['nite'] = nite
    return out

def download_votable(url):
    pass

def download_exposure(expnum,outfile=None,votable=None,dryrun=False,verbose=False):
    try:
        expnum = int(expnum)
    except ValueError:
        msg = "Exposure number could not be cast to int"
        raise ValueError(msg)

    if votable is None:
        votable = VOTABLE
        
    data = vot2npy(votable)

    match = np.where(data['expnum'] == expnum)[0]
    if len(match) == 0:
        msg = "No match to exposure: %s"%expnum
        raise Exception(msg)
    if len(match) > 1:
        msg = "Multiple matches to exposure: %s"%expnum
        raise Exception(msg)
    i = match[0]
        
    url = data[i]['reference']

    if outfile is None:
        outfile = 'DECam_{expnum:08d}.fits.fz'.format(expnum=expnum)
    cert = 'data/certificates/drlicawagnera.cert'

    cmd = noao_curl.format(url=url,outfile=outfile,cert=cert)
    if verbose: print cmd
    if dryrun: return cmd
    out = subprocess.check_output(cmd,shell=True)
    return outfile

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    print get_noao_query()
    
    for exposure in []:
        download_exposure

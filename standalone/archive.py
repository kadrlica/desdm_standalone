#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os.path
import logging
import subprocess

from standalone import noao
from standalone.templates import noao_curl, desdm_wget

DATADIR = '/data/%(drive)s/data/DTS/src'
ARCH = 'des51.b'
ARCHDIR = DATADIR%dict(drive='des51.b')
DIRNAME = os.path.join(ARCHDIR,'%(nite)8i')
FILENAME = 'DECam_%(expnum)08d.fits'

FNAL_ARCH = '/data/%(drive)s/data/DTS/src'
DESDM_ARCH = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/DTS/raw"

FILEDIR = "%(nite)8i/DECam_%(expnum)08d.fits.fz"
FILENAME = os.path.join(FNAL_ARCH,FILEDIR)

VOTABLE = '/home/s1/kadrlica/projects/standalone/data/votables/votable_20160131.vot'
NOAO_CERT = '/home/s1/kadrlica/projects/standalone/data/certificates/drlicawagnera.cert'

def get_dts_content(arch='/data/des51.b/data/src'):
    data = np.array(glob.glob(arch+'/*/DECam_*.fits.fz'))
    data = np.sort(data)
    nite,sep,exp = np.char.partition(data,'/').T
    nite = nite.astype(int)
    expnum = np.char.partition(np.char.partition(exp,'_')[:,-1],'.')[:,0].astype(int)
    return data, nite, expnum

def fill_archive(votable=None):
    if votable is None:
        votable = VOTABLE
    data = noao.vot2npy(votable)

    for d in data:
        params = dict(nite=d['nite'],expnum=d['expnum'],drive=ARCH)
        outfile = FILENAME%params
        if os.path.exists(outfile):
            msg = "Found exposure: %s"%outfile
            logging.info(msg)
            continue

        outdir = os.path.dirname(outfile)
        if not os.path.exists(outdir):
            try: os.makedirs(outdir)
            except OSError: pass

        # First try to link the file if it is already on disk
        for drive in ['des30.b','des40.b']:
            params['drive'] = drive
            archfile = FILENAME%params
            if os.path.exists(archfile):
                cmd = 'ln -s %s %s'%(archfile,outfile)
                logging.info(cmd)
                try: subprocess.check_call(cmd,shell=True)
                except subprocess.CalledProcessError: pass
                break
        if os.path.exists(outfile): continue

        # Then try to grab from DESDM
        desdm = False
        url = os.path.join(DESDM_ARCH,FILEDIR)%params
        cmd = desdm_wget.format(url=url,outfile=outfile)
        logging.info(cmd)
        try: subprocess.check_call(cmd,shell=True)
        except subprocess.CalledProcessError: pass
        if os.path.exists(outfile): continue

        # Finally try from NOAO
        cmd = 'csub '+noao_curl.format(url=d['reference'],outfile=outfile,cert=NOAO_CERT)
        logging.info(cmd)
        try: subprocess.check_call(cmd,shell=True)
        except subprocess.CalledProcessError: pass
        if os.path.exists(outfile): continue

        break
        

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    fill_archive()

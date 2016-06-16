#!/usr/bin/env python
"""
Run the DESDM pipeline steps.
"""

from configparser import SafeConfigParser, NoOptionError
import subprocess
from os.path import join

#from despyfits.DESImage import DESImage
import fitsio

import desdmlib2 as desdmlib
import dcache

def set_default(config,key,value):
    config.defaults()[key] = str(value)

def set_ccdnum(config,ccdnum):
    set_default(config,'ccdnum','%02d'%int(ccdnum))

def get_ccdnums(config):
    return ['%02d'%int(c) for c in config.defaults()['ccdnums'].split(',')] 

def print_ccdnum(ccdnum,fill='=',width=50):
    print '{:{fill}^{width}}'.format(' CCD %s '%ccdnum,fill=fill,width=width)


def copy_from_dcache(config,force=False,verbose=False):

    if isinstance(config,basestring):
        configfile = config
        config = SafeConfigParser()
        config.read(configfile)

    # Dcache Directories
    datadir = config.get('dcache','datdir')
    calbdir = config.get('dcache','caldir')
    confdir = config.get('dcache','cfgdir')
    biasdir = config.get('dcache','biasdir')
    flatdir = config.get('dcache','flatdir')
    bpmdir  = config.get('dcache','bpmdir')
    lindir  = config.get('dcache','lindir')
    bfdir   = config.get('dcache','bfdir')
    stardir = config.get('dcache','stardir')
    skydir  = config.get('dcache','skydir')

    # Data Files
    datafiles = [
        join(datadir,config.get('crosstalk','infile')),
        ]
    dcache.copy(datafiles,verbose=verbose)

    # Configuration Files
    conffiles = [
        join(confdir,config.get('crosstalk','replace')),
        join(confdir,config.get('crosstalk','crosstalk')),
        join(confdir,config.get('scamp','config')),
        join(confdir,config.get('scamp','aheader_global')),
        join(confdir,config.get('psfex','config')),
        join(confdir,config.get('updatewcs','hdupcfg')),
        ]

    # This gets a bit hairy, probably better to try to grab all 5 for
    # each sextractor run and check with dcache_ls whether the file
    # actually exists
    options = ['config','parameters_name','filter_name','starnnw_name','psf_name']
    for option in options:
        conffiles += [join(confdir,config.get('sextractor',option))]
    for option in options[:2]:
        conffiles += [join(confdir,config.get('sextractorsky',option))]
    for option in options[:4]:
        conffiles += [join(confdir,config.get('sextractorpsf',option))]

    dcache.copy(conffiles,verbose=verbose)

    # Exposure Calibration files
    calbfiles = [
        join(lindir,config.get('pixcorrect_im','lincor')),
        join(bfdir,config.get('pixcorrect_im','bf')),
        join(skydir,config.get('skyfit','pcfilename')),
        ]
    dcache.copy(calbfiles,verbose=verbose)
    
    # CCD Calibration Files
    ccdfiles = []
    for ccdnum in get_ccdnums(config):
        set_ccdnum(config,ccdnum)
        ccdfiles += [
            join(biasdir,config.get('pixcorrect_im','bias')),
            join(flatdir,config.get('pixcorrect_im','flat')),
            join(bpmdir,config.get('pixcorrect_im','bpm')),
            join(stardir,config.get('starflat','starflat')),
            join(skydir,config.get('skysubtract','pcfilename')),
            ]
    dcache.copy(ccdfiles,verbose=verbose)

    return True

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('config',help='Configuration file')
    parser.add_argument('--dryrun',action='store_true',help="Don't do anything")
    parser.add_argument('-v','--verbose',action='store_true',help="Chatter")
    args = parser.parse_args()

    dryrun = args.dryrun
    verbose = args.verbose

    config = SafeConfigParser()
    config.read(args.config)
    set_default(config,'chatter',int(verbose))

    ccdnums = get_ccdnums(config)

    ### Download the necessary data and calibrations
    copy_from_dcache(config,verbose=verbose)


    ### First run cross talk
    crosstalk = desdmlib.crosstalk()
    crosstalk.parse_config(config)
    crosstalk.run(verbose=verbose,dryrun=dryrun)

    ### Run first pixcorrect step on each CCD
    pixcorrect = desdmlib.pixcorrect()
    nullweight = desdmlib.nullweight() 
    sextractor = desdmlib.sextractor() 
    psfex      = desdmlib.psfex() 

    for ccdnum in ccdnums:
        print_ccdnum(ccdnum)
        set_ccdnum(config,ccdnum)

        pixcorrect.parse_config(config)
        pixcorrect.run(dryrun=dryrun,verbose=verbose)

        nullweight.parse_config(config)
        nullweight.run(dryrun=dryrun,verbose=verbose)

        sextractor.parse_config(config)
        sextractor.run(dryrun=dryrun,verbose=verbose)

        psfex.parse_config(config)
        psfex.run(dryrun=dryrun,verbose=verbose)

    # Re-combine catalogs
    catcombine = desdmlib.catcombine()
    catcombine.parse_config(config)
    catcombine.run(dryrun=dryrun,verbose=verbose)

    # Calculate astrometry with scamp
    scamp = desdmlib.scamp()
    scamp.parse_config(config)
    scamp.run(dryrun=dryrun,verbose=verbose)

    # Perform bleed masking and 
    updatewcs = desdmlib.updatewcs()
    bleedmask = desdmlib.bleedmask()
    skycompress = desdmlib.skycompress()
    for ccdnum in ccdnums:
        print_ccdnum(ccdnum)
        set_ccdnum(config,ccdnum)

        updatewcs.parse_config(config)
        updatewcs.run(dryrun=dryrun,verbose=verbose)

        bleedmask.parse_config(config)
        bleedmask.run(verbose=verbose,dryrun=dryrun)

        skycompress.parse_config(config)
        skycompress.run(verbose=verbose,dryrun=dryrun)

    # Create an input list
    subprocess.call('ls *bleedmask-mini.fits > listpcain',shell=True)

    skycombine = desdmlib.skycombine()
    skycombine.parse_config(config)
    skycombine.run(verbose=verbose,dryrun=dryrun)

    skyfit = desdmlib.skyfit()
    skyfit.parse_config(config)
    skyfit.run(verbose=verbose,dryrun=dryrun)

    # Run the second pass of pixcorrect
    skysubtract = desdmlib.skysubtract()

    sextractorsky = desdmlib.sextractor()
    sextractorsky.section = 'sextractorsky'

    pixcorrect_sf = desdmlib.pixcorrect()
    pixcorrect_sf.section = 'starflat'

    nullweightbkg = desdmlib.nullweight()
    nullweightbkg.section = 'nullweightbkg'

    immask = desdmlib.immask()

    rowinterp = desdmlib.rowinterp()

    sextractorpsf = desdmlib.sextractor()
    sextractorpsf.section = 'sextractorpsf'

    for ccdnum in ccdnums:
        print_ccdnum(ccdnum)
        set_ccdnum(config,ccdnum)
         
        skysubtract.parse_config(config)
        skysubtract.run(verbose=verbose,dryrun=dryrun)
         
        pixcorrect_sf.parse_config(config)
        pixcorrect_sf.run(verbose=verbose,dryrun=dryrun)
         
        nullweightbkg.parse_config(config)
        nullweightbkg.run(verbose=verbose,dryrun=dryrun)
         
        sextractorsky.parse_config(config)
        sextractorsky.run(verbose=verbose,dryrun=dryrun)
         
        immask.parse_config(config)
        immask.run(verbose=verbose,dryrun=dryrun)

        rowinterp.parse_config(config)
        rowinterp.run(verbose=verbose,dryrun=dryrun)

        # Can this be moved to the beginning of the loop?
        h = fitsio.read_header(config.get('rowinterp','out'))
        fwhm = '%0.4f'%(0.2626*float(h['FWHM']))
        config.set('sextractorpsf','seeing_fwhm',fwhm)

        sextractorpsf.parse_config(config)
        sextractorpsf.run(verbose=verbose,dryrun=dryrun)
         
        infile = config.get('sextractorpsf','catalog_name')
        outfile = infile.splitext()[0] + '.reg'
        read_geometry(infile,outfile,ccdnum)
        break

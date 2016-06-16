#!/usr/bin/env python
"""
Create the job file for all exposures in the list file.

The list file should be a csv file containing:

#nite,expnum,filter
20121124,155248,r
20121124,155249,i
...

"""

import numpy as np
from numpy import array, savetxt, recfromcsv
import subprocess
import os,shutil
from os.path import basename,dirname,join

try: 
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser
from StringIO import StringIO

import dcache

script = """#!/usr/bin/env bash
OLDHOME=$HOME
export HOME=$PWD

# Setup software from cvmfs
ulimit -a
source /cvmfs/des.opensciencegrid.org/eeups/startupcache.sh
setup fitscombine yanny
setup finalcut Y2A1dev+6
setup scamp 2.1.1+5

# Copy the code from dcache
#ifdh cp -D /pnfs/des/persistent/desdm/code/desdmLiby1e2.py .
#ifdh cp -D /pnfs/des/persistent/desdm/code/run_desdmy1e2.py .
ifdh cp -D {codedir}/{tarball} .
tar -xf {tarball}

# Write the configuration
rm -f confFile
cat <<EOF >> confFile

{config}

EOF

# Run the job
python {desdmrun} confFile >& {logfile}
du -skh .

# Copy the logfile
ifdh cp -D {logfile} {outdir}

# Remove the local file
{cleanup} 

export HOME=$OLDHOME
"""

submit = """#!/usr/bin/env bash

# Use the system python
unset PYTHONPATH
export PATH=/usr/bin/:$PATH

jobsub_submit \\
    -G des -M --disk=85GB --OS=SL6 \\
    --resource-provides=usage_model=DEDICATED \\
    file://{outbase};

"""

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('inlist',metavar='inlist.csv',
                    help="Input csv file containing: nite,expnum,filter")
parser.add_argument('-c','--config', default = None,
                    help = "Input onfigureation file")
parser.add_argument('--dryrun',action='store_true',
                    help="Do not run the job")
parser.add_argument('-r','--reqnum',default=0, type=int,
                    help="Procesing request number")
parser.add_argument('-p','--attnum',default=1, type=int,
                    help="Procesing attempt number")
parser.add_argument('-f','--force',action='store_true',
                    help="Force overwrite of output")
parser.add_argument('-q','--queue',default='local',
                    choices=['local','grid','csub'],
                    help="Queue to submit to")
args = parser.parse_args()

basedir = os.getcwd()

# Load exposures (explicitly define columns)
data = np.recfromcsv(args.inlist,names=['nite','expnum','filter'])
data = np.atleast_1d(data)

local = True if args.queue=='local' else False
params = dict(reqnum=args.reqnum,attnum=args.attnum)

config = ConfigParser()
if args.config is None:
    codedir = os.path.dirname(os.path.abspath(__file__))
    configfile = os.path.join(codedir,'config','default.cfg')
else:
    configfile = args.config
config.read(configfile)
 
# Set reqnum and attnum
config.set('General','r','{reqnum:04d}'.format(**params))
config.set('General','p','{attnum:02d}'.format(**params))
config.set('general','reqnum','{reqnum:04d}'.format(**params))
config.set('general','attnum','{attnum:02d}'.format(**params))

# Copy code to dcache
desdmrun = config.get('general','desdmrun')
version = config.get('general','version')
srcdir = config.get('general','srcdir').rstrip('/')
codedir = config.get('dcache','srcdir').rstrip('/')

# Create a tarball of the current version of the code
os.chdir(srcdir)
tarball ='standalone-%s.tar'%version
cmd = "git archive %s -o %s"%(version,join(basedir,tarball))
print cmd
subprocess.check_output(cmd,shell=True)
os.chdir(basedir)

# Make the code directory in dcache and copy the tarball
dcache.mkdir(codedir,'-p',verbose=True)
dcache.cp(join(basedir,tarball),codedir,)

params.update(codedir=codedir,desdmrun=desdmrun,tarball=tarball)

print "Preparing %s exposures..."%len(data)
for i,d in enumerate(data):
    print 30*"-"
    os.chdir(basedir)

    params.update(zip(d.dtype.names,d))

    # create the exposure directory 
    rundir = '{expnum}'.format(**params)
    if os.path.exists(rundir):
        if not args.force:
            print "Found %s; skipping..."%rundir
            continue
        else:
            shutil.rmtree(rundir)

    os.makedirs(rundir)
    
    config.set('General','nite',params['nite'])
    config.set('General','expnum',params['expnum'])
    config.set('General','filter',params['filter'])
    config.set('general','nite',params['nite'])
    config.set('general','expnum',params['expnum'])
    config.set('general','filter',params['filter'])
    
    configio = StringIO()
    config.write(configio)
    params['config'] = configio.getvalue()

    # Add some more 
    params['logfile'] = "out{expnum}.log".format(**params)
    params['outdir'] = "/pnfs/des/scratch/kadrlica".format(**params) #Testing

    cleanup=('#' if local else '')
    cleanup+="rm -f *.fits *.fits.fz *.ps *.psf *.xml full_1.cat *.head"
    params['cleanup'] = cleanup

    copylog=('#' if local else '')
    copylog+="ifdh cp -D {logfile} {outdir}".format(**params)
    params['copylog'] = copylog

    # Write the run script
    outbase = 'job_D{expnum:08d}_r{reqnum:d}p{attnum:02d}.sh'.format(**params)
    outfile = os.path.join(rundir,outbase)
    out = open(outfile,'w')
    out.write(script.format(**params))
    out.close()
    os.chmod(outfile,0755)

    # Write the submission script
    subfile = os.path.join(rundir,'submit.sh')
    out=open(subfile,'w')
    out.write(submit.format(outbase=outbase))
    out.close()
    os.chmod(subfile,0755)

    # Run the script
    if args.queue == 'local':
        # Run the script right here...
        cmd = "sh %s"%outbase
    elif args.queue == 'csub':
        # Run the script as a process on the local machine
        cmd = "csub sh %s"%outbase
    elif args.queue == 'grid':
        # Submit to the grid
        cmd = "sh %s"%submit

    if args.dryrun: continue
    os.chdir(rundir)
    subprocess.call(cmd,shell=True)

#!/usr/bin/env python
"""
Module with templates for various scripts
"""
__author__ = "Alex Drlica-Wagner"

unitname = "{unitname:>08s}"
reqnum   = "{reqnum:>04s}"
attnum   = "{attnum:>02s}"

desdmjob = """#!/usr/bin/env bash
OLDHOME=$HOME
export HOME=$PWD

# Setup the software
ulimit -a
source /cvmfs/des.opensciencegrid.org/eeups/startupcache.sh
setup fitscombine yanny
setup finalcut Y2A1dev+6
setup scamp 2.1.1+5

#ifdh cp -D /pnfs/des/persistent/desdm/code/desdmLiby1e2.py .
#ifdh cp -D /pnfs/des/persistent/desdm/code/run_desdmy1e2.py .
ifdh cp -D {codedir}/{desdmlib} .
ifdh cp -D {codedir}/{desdmrun} .

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

### JobSub ###
jobsubmit = """#!/usr/bin/env bash

# Use the system python
unset PYTHONPATH
export PATH=/usr/bin/:$PATH

jobsub_submit \\
    -G des -M --disk=85GB --OS=SL6 \\
    --resource-provides=usage_model=DEDICATED \\
    file://{outbase};

"""

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

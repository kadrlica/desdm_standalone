#!/usr/bin/env python
"""
Python wrapper to work with Fermilab dcache.

Note: There is a SWIG wrapped pyhton interface to ifdh:
/cvmfs/fermilab.opensciencegrid.org/products/common/prd/ifdhc/v1_8_2/Linux64bit-2-6-2-12/lib/python/ifdh.py

"""
__author__ = "Alex Drlica-Wagner"
__all__ = ['dcache_cp','dcache_mkdir','dcache_rm','dcache_ls','dcache_chmod',
           'decache_mv']

import os
import subprocess

def parse_opts(opts):
    if opts is None:
        return []
    elif isinstance(opts,basestring):
        return opts.split()
    else:
        return opts
    msg = 'Unrecognized options: %s'%opts
    raise Exception(msg)

def execute(cmd,dryrun=False,verbose=False):
    if verbose: print cmd
    if dryrun: return cmd
    return subprocess.check_output(cmd,shell=True)

def dcache_isdir(path):
    pathbase = os.path.basename(path)

    if os.path.isdir(path): 
        return True
    if os.path.isfile(path):
        return False
    if pathbase in ['','.','..']: 
        return True

    # This is more expensive, so do it last...
    try: 
        out = dcache_ls(path).strip().split('\n')
        # Easy if empty dir or contains more than one file
        if len(out) == 0 or len(out) > 1:
            return True

        outbase = os.path.basename(out[0])
         # This is harder...
        if outbase == '':
            # For empty directories, idfh returns the trailing '/'
            return True
        if outbase != pathbase:
            # The directory and it's contents have different names
            return True

    except subprocess.CalledProcessError:
        # Path doesn't exist (new file)
        pass

    return False

def dcache_copy(src,dest=None,opts=None,dryrun=False,verbose=False):
    opts = parse_opts(opts)

    # If destination not specified assume current directory
    if dest is None: dest = '.'

    # Allow a list of files to copy
    if not isinstance(src,basestring):
        src = ' '.join(src)

    # Try to determine if 'dest' is a directory
    if '-D' not in opts and (len(src.split())>1 or dcache_isdir(dest)): 
        opts += ['-D']

    if '-D' in opts and '-r' in opts:
        return dcache_cpdr(src,dest,dryrun,verbose=verbose)
    elif '-D' in opts:
        return dcache_cpd(src,dest,dryrun,verbose=verbose)
    elif '-r' in opts:
        return dcache_cpr(src,dest,dryrun,verbose=verbose)

    return dcache_cp(src,dest,dryrun,verbose=verbose)

def dcache_cp(src,dest,dryrun=False,verbose=False):
    cmd = 'ifdh cp %s %s'%(src,dest)
    return execute(cmd,dryrun,verbose)

def dcache_cpr(src,dest,dryrun=False,verbose=False):
    cmd = 'ifdh cp -r %s %s'%(src,dest)
    return execute(cmd,dryrun,verbose)

def dcache_cpd(src,dest,dryrun=False,verbose=False):
    cmd = 'ifdh cp -D %s %s'%(src,dest)
    return execute(cmd,dryrun,verbose)

def dcache_cpdr(src,dest,dryrun=False,verbose=False):
    cmd = 'ifdh cp -D -r %s %s'%(src,dest)
    return execute(cmd,dryrun,verbose)

def dcache_mkdir(dirname,opts=None,dryrun=False,verbose=False):
    opts = parse_opts(opts)
    if '-p' in opts:
        return dcache_mkdir_p(dirname,dryrun)

    cmd = 'ifdh mkdir %s'%(dirname)
    return execute(cmd,dryrun,verbose)

def dcache_mkdir_p(dirname,dryrun=False,verbose=False):
    cmd = 'ifdh mkdir_p %s'%(dirname)
    return execute(cmd,dryrun,verbose)

def dcache_rm(filename,opts=None,dryrun=False,verbose=False):
    opts = parse_opts(opts)
    if '-r' in opts:
        maxdepth = 1000 # Probably large enough...
        filelist = dcache_ls('%s %i'%(filename,maxdepth)).strip().split('\n')
        fileparts = sorted([ f.split('/') for f in filelist ])[::-1]
        out = []
        for f in fileparts:
            o = dcache_rm('/'.join(f),dryrun=dryrun,verbose=verbose)
            if o: out.append(o)
        return '; '.join(out)

    cmd = 'ifdh rm %s'%filename
    return execute(cmd,dryrun,verbose)

def dcache_rmdir(dirname,dryrun=False,verbose=False):
    cmd = 'ifdh rmdir %s'%dirname
    return execute(cmd,dryrun,verbose)

def dcache_ls(filename,opts=None,dryrun=False,verbose=False):
    opts = parse_opts(opts)
    if '-s' in opts:
        return dcache_lss(filename,dryrun=dryrun,verbose=verbose)
    elif '-l' in opts:
        return dcache_ll(filename,dryrun=dryrun,verbose=verbose)

    cmd = 'ifdh ls %s'%filename
    # Should sort the output
    return execute(cmd,dryrun,verbose)

def dcache_ll(filename,dryrun=False,verbose=False):
    cmd = 'ifdh ll %s'%filename
    return execute(cmd,dryrun,verbose)

def dcache_lss(filename,dryrun=False,verbose=False):
    cmd = 'ifdh lss %s'%filename
    return execute(cmd,dryrun,verbose)

def dcache_chmod(filename,mode,dryrun=False,verbose=False):
    cmd = 'ifdh chmod %s %s'%(mode,filename)
    return execute(cmd,dryrun,verbose)

def dcache_mv(src,dest,dryrun=False,verbose=False):
    cmd = 'ifdh mv %s %s'%(src,dest)
    return execute(cmd,dryrun,verbose)

def copy_from_dcache(filename):
    return dcache_cpd(filename,'.',verbose=True)

def copy_to_dcache(filename, dirname):
    return dcache_cpd(filename,dirname,verbose=True)
    
# Some quick aliases
rm = dcache_rm
copy = dcache_copy
cp = dcache_cp
cpr = dcache_cpr
cpd = dcache_cpd
mkdir = dcache_mkdir
ls = dcache_ls
chmod = dcache_chmod
isdir = dcache_isdir

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    ##############################
    #print "Setup paths..."
    basedir = '/pnfs/fnal.gov/usr/des/scratch/kadrlica/'
    subdir = os.path.join(basedir,'hello')
    outdir = os.path.join(subdir,'world')
    # A real bastard for testing
    devilfile = os.path.join(outdir,os.path.basename(outdir))

    filename = 'tmp.txt'
    f = open(filename,'w')
    f.write('hello world')
    f.close()
    filename2 = filename+'2'

    ##############################
    print "### Checking ls... ###"
    for opt in ['','-l','-s']:
        print dcache_ls(basedir,opt,verbose=True)

    ##############################
    print "### Checking mkdir... ### "
    print dcache_mkdir(outdir,'-p',verbose=True)
    print dcache_ls(outdir,verbose=True)
    
    print dcache_mkdir(outdir,verbose=True)
    print dcache_ls(outdir,verbose=True)

    ##############################
    print "### Checking cp... ###"
    print dcache_cp(filename,outdir,verbose=True)
    print dcache_cp(os.path.join(outdir,filename),filename2,verbose=True)
    print dcache_cp(filename2,outdir,verbose=True)
    print dcache_cp([filename,filename2],outdir,verbose=True)
    print dcache_cp('.',[filename,filename2],verbose=True)

    ##############################
    print "### Checking isdir... ###"
   
    for f in [basedir,outdir,os.path.join(outdir,filename),os.getcwd(),filename]:
        print '%s :  %s'%(f,dcache_isdir(f))

    # Now this is a real doosy...
    dcache_cp(filename,devilfile,verbose=True)
    print 'DEVIL FILE :  %s'%(dcache_isdir(devilfile))
    dcache_rm(os.path.join(outdir,filename),verbose=True)
    print

    ##############################
    print "### Checking rm... ###"

    print dcache_ls(subdir+' 1000',verbose=True)
    print dcache_rm(subdir,'-r',verbose=True)
    print dcache_ls(os.path.dirname(subdir),verbose=True)

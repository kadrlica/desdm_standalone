#!/usr/bin/env python
"""
Processing steps for star flats.  Assume that headers have already been updated.
Will need to have already:
setup imdetrend
setup pixcorrect
setup scamp
setup sextractor
setup psfex
setup python
setup fitsio
"""

import sys
import os
import subprocess
import shutil
import fitsio
import socket
import re
from pixcorrect import decaminfo


# Configuration files:
config_dir = "/data3/DECAM/PHOTCAL/CONFIG"
sexstar_config = os.path.join(config_dir,"sexstar.config")
sexstar_param = os.path.join(config_dir,"sexstar.param")
sex_conv = os.path.join(config_dir,"gauss_3.0_7x7.conv")
sex_nnw = os.path.join(config_dir,"sex.nnw")
sexfinal_config = os.path.join(config_dir,"sexstarflat.config")
sexfinal_param = os.path.join(config_dir,"sexfinal.param")
scamp_config  = os.path.join(config_dir,"scamp.config")
psfex_config  = os.path.join(config_dir,"psfex.config")

# Calibration files:
cal_dir = "/data3/DECAM/PHOTCAL/CALS"
xtalk_file = os.path.join(cal_dir, "DECam_20130606.xtalk")
hupdate_file = os.path.join(cal_dir, "20140709_DES_header_update.20140303")
linearity_file = os.path.join(cal_dir,'linearity_table_v0.5.fits')
bpm_file = os.path.join(cal_dir,'20140901t0928_Y2T2/20140901t0928_Y2T2_bpm_{ccd:02d}.fits')
bf_file = os.path.join(cal_dir,'bfmodel_20150305.fits')
bias_file = os.path.join(cal_dir,'D_n20141204t1209_c{ccd:02d}_r1426p08_biascor.fits')
dome_file = os.path.join(cal_dir,'D_n20141204t1209_{band:s}_c{ccd:02d}_r1426p08_norm-dflatcor.fits')
sky_template = os.path.join(cal_dir,'sky_{band:s}_{ccd:02d}.fits')
sky_pca = os.path.join(cal_dir,'pca_{band:s}.fits')

# Files for intermediate products
raw_file = 'DECam_{expid:>08d}.fits.fz'         #Raw data
log_file = 'D{expid:08d}_{ccd:02d}.log'         #processing logs

xt_template = "D{expid:>08d}_%02d.xt.fits"      #Crosstalk-corrected
xt_file = 'D{expid:08d}_{ccd:02d}.xt.fits'
fc_file = 'D{expid:08d}_{ccd:02d}.fc.fits'      #Flat-corrected
bm_file = 'D{expid:08d}_{ccd:02d}.bm.fits'      #Bleedmasked
cr_file = 'D{expid:08d}_{ccd:02d}.cr.fits'      #CR-detected 
ss_file = 'D{expid:08d}_{ccd:02d}.ss.fits'      #sky subtracted
tmp_file = 'D{expid:08d}_{ccd:02d}.tmp.fits'    #temporary sextractor input image

mini_file = 'D{expid:08d}_{ccd:02d}.mini.fits'  #single-ccd mini-sky
mini_template = 'D{expid:08d}_%%02d.mini.fits'  

sky_file = 'D{expid:08d}.mini.fits'             #combined mini-sky
residsky_file = 'D{expid:08d}.resid.fits'       #residual sky

starcat_file = 'D{expid:08d}_{ccd:02d}.starcat' #Star-finding catalog
# This name is generated by PSFEx based on the input catalog name:
psf_file = 'D{expid:08d}_{ccd:02d}.psf'         #PSFEx output
merge_file = 'D{expid:08d}.mergecat.fits'       #SCAMP input merged catalog
head_file = 'D{expid:08d}.mergecat.head'        #SCAMP output head file
cat_file = 'D{expid:08d}_{ccd:02d}.cat'         #single-ccd output catalogs
final_file = 'D{expid:08d}.cat'                 #final merged catalog

def getScratchDir(unique):
    # Return path to a suitable scratch directory on this node
    hostname = socket.gethostname()
    match = re.search(r'(\d+)$',hostname)
    if not match:
        print "Could not find node number"
        sys.exit(1)
    if int(match.group()) == 19:
        # if False:
        # No working scratch dir:
        return os.path.join(os.getcwd(), 'tmp.' + unique)
    else:
        return '/scratch/garyb/tmp.' + unique

def detrend(expid, band, **kwargs):
    """
    From raw image, perform xtalk and pixcorrect operations through flat-field.
    """
    d = {'expid':expid, 'band':band}

    cmd = 'DECam_crosstalk ' + raw_file.format(**d) + \
          ' ' + xt_template.format(**d) +\
          ' -crosstalk ' + xtalk_file + \
          ' -replace ' + hupdate_file + \
          ' -overscanfunction 0 -overscansample 1 -overscantrim 5 ' + \
          ' -photflag 1 -verbose 0'

    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)


    for ccd in range(1,63):
        d['ccd']=ccd
        if ccd not in (2,61):
            cmd = 'pixcorrect_im -i ' + xt_file.format(**d) + \
                  ' -o ' + fc_file.format(**d) + \
                  ' -l ' + log_file.format(**d) + \
                  ' --bias ' + bias_file.format(**d) + \
                  ' --lincor ' + linearity_file + \
                  ' --gain ' + \
                  ' --bpm ' + bpm_file.format(**d) + \
                  ' --bf ' + bf_file + \
                  ' --flat ' + dome_file.format(**d) + \
                  ' --mini ' + mini_file.format(**d) + \
                  ' --resaturate'
            retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
            if retval != 0:
                sys.exit(1)
            os.remove(xt_file.format(**d))
 
    # Combine the mini-sky images
    explog_file = 'D{expid:08d}.log'
    cmd = 'sky_combine --miniskyfiles ' + mini_template.format(**d) + \
          ' --outfilename ' + sky_file.format(**d) + \
          ' -l ' + explog_file.format(**d)
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        exit(1)
    for ccd in range(1,63):
        if ccd not in (2,61):
            d['ccd']=ccd
            os.remove(mini_file.format(**d))
    return

def skysubtract(expid, band, in_file=bm_file, **kwargs):
    """
    Fit and subtract the sky, and create weight plane
    """
    d = {'expid':expid, 'band':band}

    explog_file = 'D{expid:08d}.log'
    cmd = 'sky_fit --infilename ' + sky_file.format(**d) + \
          ' --outfilename ' + residsky_file.format(**d) + \
          ' --pcfilename ' + sky_pca.format(**d) + \
          ' -l ' + explog_file.format(**d)
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        exit(1)

    for ccd in range(1,63):
        d['ccd']=ccd
        if ccd not in (2,61):
            cmd = 'sky_subtract -i ' + in_file.format(**d) + \
                  ' -o ' + ss_file.format(**d) + \
                  ' --pcfilename ' + sky_template.format(**d) + \
                  ' --fitfilename ' + residsky_file.format(**d) + \
                  ' --domefilename ' + dome_file.format(**d) + \
                  ' --weight sky ' + \
                  ' -l ' + log_file.format(**d)
            retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
            if retval != 0:
                sys.exit(1)
            os.remove(in_file.format(**d))

def bleed(expid, ccd, in_file=fc_file,**kwargs):
    """
    Run mkbleedmask. 
    """
    d = {'expid':expid, 'ccd':ccd}
    cmd = 'mkbleedmask ' + in_file.format(**d) + \
          ' ' + bm_file.format(**d) + \
          ' -m -b 5 -f 1.0 -l 7 -n 7 -r 5 -s 40 -t 20 -v 3 -w 1.5 -y 1.0 '
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)
    os.remove(in_file.format(**d))
    return

def cr(expid, ccd, in_file=bm_file):
    """
    Run cosmic ray / streak masking
    """
    d = {'expid':expid, 'ccd':ccd}
    cmd = 'immask ' + in_file.format(**d) + \
          ' ' + cr_file.format(**d) + \
          ' ???? need to set this all up!!!'
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)
    os.remove(in_file.format(**d))
    return

def sexstar(expid, ccd, band, in_file=ss_file, **kwargs):
    """
    Run SExtractor on single image to extract catalog for scamp and psfex
    """
    d = {'expid':expid, 'ccd':ccd, 'band':band}

    # Preprocess the image for SExtractor by nulling weights in bad regions
    cmd = 'null_weights -i ' + in_file.format(**d) + \
                  ' -o ' + tmp_file.format(**d) + \
                  ' -l ' + log_file.format(**d) + \
                  ' --bitmask BADAMP,EDGEBLEED,STREAK,EDGE' + \
                  ' --resaturate'
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)

    cmd = 'sex ' + tmp_file.format(**d)+'[SCI]' + \
          ' -c '  + sexstar_config + \
          ' -PARAMETERS_NAME ' + sexstar_param + \
          ' -FILTER_NAME ' + sex_conv + \
          ' -STARNNW_NAME ' + sex_nnw + \
          ' -CATALOG_NAME ' + starcat_file.format(**d) + \
          ' -FLAG_IMAGE ' +  tmp_file.format(**d)+'[1]' + \
          ' -WEIGHT_IMAGE ' +  tmp_file.format(**d)+'[2]' + \
          ' -BACKPHOTO_TYPE LOCAL'

    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)
    os.remove(tmp_file.format(**d))
    return
    
def mergeLDAC(single_template, outfile, info, cleanup=False):
    """
    Concatenate LDAC-format catalogs from single ccds into a single catalog with all
    CCDs.
    single_template is a string that can be formatted with the fields in the
    info  dictionary, augmented with a 'ccd' key, to form names of individual catalogs.
    outfile is the name of the merged catalog.

    cleanup=True deletes constituent catalogs

    Returns a list of the CCDNUMs that were found and concatenated.
    """

    mergecat = None

    ccds_used = []
    for ccd in range(1,63):
        info['ccd'] = ccd
        if os.path.exists(single_template.format(**info)):
            ccds_used.append(ccd)
            if mergecat is None:
                # Copy the first file to become the merged catalog
                shutil.copyfile(single_template.format(**info), outfile)
                mergecat = fitsio.FITS(outfile, 'rw')
            else:
                # Append 2 extensions to the merged catalog
                for extn in ('LDAC_IMHEAD','LDAC_OBJECTS'):
                    data,header = fitsio.read(single_template.format(**info),ext=extn,header=True)
                    mergecat.write(data, header=header,extname=extn)
    mergecat.close()
    if cleanup:
        for ccd in range(1,63):
            info['ccd'] = ccd
            if os.path.exists(single_template.format(**info)):
                os.remove(single_template.format(**info))
        
    return ccds_used

def scamp(expid, in_file=ss_file, **kwargs):
    """
    Concatenate catalogs from an exposure, run scamp on them vs
    reference catalog.
    (then put solution back into image headers)
    """
    # First merge the LDAC catalogs
    d = {'expid':expid}
    ccds_used = mergeLDAC(starcat_file, merge_file.format(**d),d,
                          cleanup=False)
    
    # run scamp
    cmd = 'scamp ' + merge_file.format(**d) + \
          ' -c ' + scamp_config
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)
        
    # Suck WCS back into image headers
    skip = ('COMMENT','HISTORY','')
    ccds_used.reverse()
    d['ccd'] = ccds_used.pop()
    hdr = fitsio.FITSHDR()
    heads = open(head_file.format(**d))
    for line in heads:
        if line.split()[0].strip()=='END':
            # We've reached the end of one image's header.  Transfer values to this
            # ccd's image header and move to next CCD
            with fitsio.FITS(in_file.format(**d),'rw') as fits:
                for k in hdr.keys():
                    if k not in skip:
                        fits[0].write_key(k, hdr[k])
            if len(ccds_used)>0:
                d['ccd'] = ccds_used.pop()
            else:
                d['ccd'] = None
            hdr = fitsio.FITSHDR()
        else:
            hdr.add_record(line.strip())
    # Clean up
    os.remove(merge_file.format(**d))
    os.remove(head_file.format(**d))
    return

def catalog(expid, ccd, band, in_file=ss_file, cleanup=False, **kwargs):
    """
    Make a catalog, including PSF fitting to get SPREAD_MODEL
    """
    
    d = {'expid':expid, 'ccd':ccd, 'band':band}

    cmd = 'psfex ' + starcat_file.format(**d) + \
          ' -c '  + psfex_config
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)

    if cleanup:
        os.remove(starcat_file.format(**d))

    # Preprocess the image for SExtractor by nulling weights in bad regions
    cmd = 'null_weights -i ' + in_file.format(**d) + \
                  ' -o ' + tmp_file.format(**d) + \
                  ' -l ' + log_file.format(**d) + \
                  ' --bitmask BADAMP,EDGEBLEED,STREAK,EDGE' 
    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)

    cmd = 'sex ' + tmp_file.format(**d)+'[SCI]' + \
          ' -c '  + sexfinal_config + \
          ' -PARAMETERS_NAME ' + sexfinal_param + \
          ' -FILTER_NAME ' + sex_conv + \
          ' -STARNNW_NAME ' + sex_nnw + \
          ' -PSF_NAME ' + psf_file.format(**d) + \
          ' -CATALOG_NAME ' + cat_file.format(**d) + \
          ' -FLAG_IMAGE ' +  tmp_file.format(**d)+'[1]' + \
          ' -WEIGHT_IMAGE ' +  tmp_file.format(**d)+'[2]' + \
          ' -BACKPHOTO_TYPE LOCAL'

    retval = subprocess.call(cmd.split(),stderr=subprocess.STDOUT)
    if retval != 0:
        sys.exit(1)
    os.remove(psf_file.format(**d))
    os.remove(tmp_file.format(**d))
              
    return

if __name__=='__main__':
    expid = int(sys.argv[1])

    # Make a subdirectory for the processed images
    try:
        os.mkdir('IMAGES')
    except OSError:
        print "IMAGES directory couldn't be created, continuing"
    try:
        os.mkdir('CATS')
    except OSError:
        print "CATS directory couldn't be created, continuing"
    
    # Move to the working directory
    startDir = os.getcwd()
    scratchDir = getScratchDir('{:>08d}'.format(expid))
    os.makedirs(scratchDir)
    os.chdir(scratchDir)

    # Link the raw image into the scatch directory
    d = {'expid':expid}
    raw = raw_file.format(**d)
    os.symlink(os.path.join(startDir,'RAW/',raw),raw)
    # Get band name from FILTER keyword
    band = decaminfo.get_band(fitsio.read_header(raw,0)['FILTER'])
    d['band'] = band
    
    # Do detrending
    detrend(**d)
    # Do bleed/star masking 
    for ccd in range(1,63):
        if ccd not in (2,61):
            d['ccd']=ccd
            bleed(**d)

    # Substract sky, create weights
    skysubtract(**d)

    # produce star catalogs
    for ccd in range(1,63):
        if ccd not in (2,61):
            d['ccd']=ccd
            sexstar(**d)

    # Astrometric shift solution
    scamp(**d)

    # Final catalogs with PSF fitting
    for ccd in range(1,63):
        if ccd not in (2,61):
            d['ccd']=ccd
            catalog(cleanup=True, **d)
    # Combine catalogs into one
    mergeLDAC(cat_file, final_file.format(**d), d)
    
    # Copy the catalog and images back to base directory (and miniskies)
    for ccd in range(1,63):
        d['ccd']=ccd
        src = ss_file.format(**d)
        if os.path.exists(src):
            shutil.move(src, os.path.join(startDir,'IMAGES/',src))
    src = sky_file.format(**d)
    shutil.move(src, os.path.join(startDir,'IMAGES/',src))
    src = residsky_file.format(**d)
    shutil.move(src, os.path.join(startDir,'IMAGES/',src))
    src = final_file.format(**d)
    shutil.move(src, os.path.join(startDir,'CATS/',src))

    os.chdir(startDir)
    # Torch the work area:
    shutil.rmtree(scratchDir)
    sys.exit(0)

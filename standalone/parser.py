#!/usr/bin/env python
import argparse
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import subprocess

mkbleedhelp = """
Detects bleedtrails in Y for incoming FITS image.
Sets BADPIX_TRAIL for detected bleedtrails,
and BADPIX_INTERP if interpolation is enabled,
and BADPIX_STAR if star masking is enabled.


mkbleedmask Usage: 

mkbleedmask [-aeghimzD] [-b <factor> -d <npix> -f <value> -j <npix> -l <value> -n <npix> -o <filename> -r <factor> -s <n> -t <npixels> -v <level> -w <value> -x <filename> -y <value> -E <num_x_dilations> -L <npix> -W <factor> ] <infile> <outfile> 

-a,--interpolate-stars
Do radial interpolation over saturated stars for which bleedtrails were detected.

-b,--bgreject <factor>
Use specified <factor> as scalefactor for background rejection. (5.0)

-d,--edgesize <npix>
Size of edges used to detect edgebleed. (15)

-e,--version
Print version and exit

-f,--scalefactor <value>
Use specified <value> as scalefactor on background to detect bleedtrails. (1.0)

-g,--global_only
Do not use local statistics.

-h,--help
Prints this long version of help.

-i,--interpolate
Interpolate over bleedtrails.

-j,--trailreject <npix>
Reject bleedtrails of size <npix> or smaller. (1)

-l,--starlevel <value>
Use specified <value> as scalefactor on background to detect stars. (5.0)

-m,--starmask
Create a mask for the detected bright objects. (No)

-n,--numtrail <npix>
Number of contiguous pixels required for trail detection. (trail_length/2)

-o,--saturated_objects <filename>
Filename for saturated star table.(None)

-r,--growbox <factor>
Factor by which to grow blob boxes for determining background level. (2)

-s,--bgiters <n>
Number of iterations for background estimation. (1)

-t,--trail_length <npixels>
Use specified <npixels> for trail structuring element in Y. (20)

-v,--verbose <level>
Verbosity level [0-3], default = 1.

-w,--growrad <value>
Factor by which to grow detected star radii. (1.0)

-x,--trailboxes <filename>
Filename for trail box table.(None)

-y,--holelevel <value>
Number of sigma below sky for hole detection. (3.0)

-z,--zerostarweights
Set weights to zero for masked stars.

-D,--Debug
Turn on debug/development mode.

-E,--Expand <num_x_dilations>
Expand bleed trails using dilation in order to capture anomalous adjacent pixels (0).

-L,--Long_strong <npix>
Length scale used to check when dilating long/strong bleed trails (30).

-W,--trail_weight_factor <factor>
Factor by which to downweight pixels flagged as bleed trails. (0)

<infile>
Input FITS file.


<outfile>
Output FITS file.
"""

xtalkhelp = """
DECam_crosstalk <infile.fits> <outfile> <options>
  -crosstalk <crosstalk matrix- file>
  -linear 
  -crossatthresh <factor> 
  -photflag <0 or 1>
  -presatmask
  -postsatmask
  -overscan
  -overscansample <-1 for MEDIAN, 0 for MEAN, 1 for MEAN w/MINMAX>, 
  -overscanfunction < -50 < N < -1 for cubic SPLINE, 0 for LINE_BY_LINE, 1 < N < 50 for legendre polynomial>
  -overscanorder <1-6>, default=1 <order for legendre polynomial>
  -overscantrim <ncols>, default=0 <trim ncols at both edges of overscan, default = 0>
  -maxhdunum <integer>, default=62 <all valid images>
  -ccdlist <comma-separated list of CCDs to process>
  -hdulist <comma-separated list of HDUs to process>
  -focuschipsout
  -replace <filename>
  -verbose <0-3>
  -help (print help message and exit)
  -version (print version and exit)
"""

class DESDMApp(object):
    _prog = "desdm_app"
    _description = 'Python wrapper around %s.'%_prog
    _section = 'desdmapp'

    def __init__(self,options=None,**kwargs):
        self.namespace = argparse.Namespace()
        self.parser = self._create_parser()
        self.namespace = self.parser.parse_args(options,self.namespace)
        self.namespace.__dict__.update(**kwargs)

    def __call__(self,dryrun=False):
        return self.run(dryrun=dryrun)

    def _create_parser(self):
        parser = ArgumentParser(prog=self._prog,description=self._description,
                                formatter_class=ArgumentDefaultsHelpFormatter,
                                add_help=False)

        for args,kwargs in self._arguments:
            parser.add_argument(*args,**kwargs)
            
        return parser

    def run(self,dryrun=False):
        cmdline = self.cmdline()
        if dryrun: return cmdline
        return subprocess.check_output(cmdline,shell=True)

    def parse_args(self,args=None):
        namespace = getattr(self,'namespace',None)
        self.namespace = self.parser.parse_args(args,namespace)
        return self.namespace
        
    def parse_config(self,config,section=None):
        if section is None: section = self._section

        if isinstance(config,basestring):
            configfile = config
            config = SafeConfigParser()
            config.read(configfile)

        for key,value in config.items(section):
            setattr(self.namespace,key,value)
        return self.namespace

    def cmdlist(self):
        cmdlist = [self._prog]
        for action in self.parser._actions:

            dest = action.dest
            default = action.default
            value = getattr(self.namespace,dest)
            if value == default: continue
         
            if not action.option_strings: option = dest
            else: option = action.option_strings[0]

            if not action.option_strings:
                # Positional argument
                cmdlist += [value]
            elif isinstance(action,argparse._StoreTrueAction) and value:
                cmdlist += [option]
            elif isinstance(action,argparse._StoreFalseAction) and not value:
                cmdlist += [option]
            elif isinstance(action,argparse._HelpAction):
                cmdlist += [option]
            elif isinstance(action,argparse._CountAction):
                cmdlist += value*[option]
            elif isinstance(action,argparse._StoreAction):
                cmdlist += [option,value]
            else:
                msg = "Invalid action type: %s"%type(action)
                raise Exception(msg)

        return cmdlist

    def cmdline(self):
        cmdlist = self.cmdlist()
        return ' '.join([str(x) for x in cmdlist])

class crosstalk(DESDMApp):
    _prog = 'DECam_crosstalk'
    _section = 'crosstalk'

    def _create_parser(self):
        """
        Create an ArgumentParser from help message.
        """
        parser = ArgumentParser(prog=self._prog,description=self._description,
                                formatter_class=ArgumentDefaultsHelpFormatter,
                                add_help=False)

        cmd = '%s -help'%self._prog
        try:
            out = subprocess.check_output(cmd,shell=True)
        except subprocess.CalledProcessError:
            out = xtalkhelp
        lines = out.strip().split('\n')
        arguments = ['<infile>','<outfile>']
        options = [l.strip() for l in lines if l.strip().startswith('-')]

        for arg in arguments:
            args = (arg.strip('<>'),)
            kwargs = dict(nargs='?')
            parser.add_argument(*args,**kwargs)

        for opt in options:
            olist = opt.split()
            if (len(olist) == 1) or (olist[1][0] != '<'): 
                args = olist[0].split(',')
                kwargs = dict(action='store_true')
            else:
                args = (olist[0],)
                kwargs = dict(action='store')
            parser.add_argument(*args,**kwargs)

        return parser


class mkbleedmask(DESDMApp):
    _prog = 'mkbleedmask'
    _section = 'bleedmask'

    def _create_parser(self):
        """
        Create an ArgumentParser from 'mkbleedmask -h'.
        """
        parser = ArgumentParser(prog=self._prog,description=self._description,
                                formatter_class=ArgumentDefaultsHelpFormatter,
                                add_help=False)

        cmd = '%s -h'%self._prog
        try: 
            out = subprocess.check_output(cmd,shell=True)
        except subprocss.CalledProcessError:
            out = mkbleedhelp
        lines = [l.strip() for l in out.strip().split('\n')]
        options = [l for l in lines if l.startswith('-')]

        arguments = [l for l in lines if l.startswith('<') and l.endswith('>')]

        for arg in arguments:
            args = (arg.strip('<>'),)
            kwargs = dict(nargs='?')
            parser.add_argument(*args,**kwargs)

        for opt in options:
            olist = opt.split()
            if len(olist) == 1: 
                args = olist[0].split(',')
                kwargs = dict(action='store_true')
            else:
                args = olist[0].split(',')
                kwargs = dict(action='store')
            parser.add_argument(*args,**kwargs)

        return parser

if __name__ == '__main__':

    # cmdtest = 'DECam_crosstalk DECam_00229650.fits.fz D00229650_g_%02d_r0000p01_xtalk.fits -crosstalk DECam_20130606.xtalk -ccdlist 3,4 -overscanfunction 0 -overscansample 1 -overscantrim 5 -photflag 1 -verbose 0 -replace DES_header_update.20151120'
    # app = crosstalk()
    # app.parse_args(cmdtest.split()[1:])
    # print cmdtest
    # print app(dryrun=True)
    # print

    # cmdtest='pixcorrect_im --verbose --in xtalk.fits -o detrend.fits --bias biascor.fits --bpm bpm.fits --lincor lin.fits --bf bfmodel.fits --gain  --flat dflatcor.fits --resaturate --fixcols --addweight'
    # app = pixcorrect()
    # app.parse_args(cmdtest.split()[1:])
    # print cmdtest
    # print app(dryrun=True)
    # print
    
    # cmdtest='mkbleedmask inname.fits outname.fits -m -b 5 -f 1.0 -l 7 -n 7 -r 5 -t 20 -v 3 -w 2.0 -y 1.0 -s 100 -v 3 -E 6 -L 30 -x trailbox.fits -o satstars.fits'
    # app = mkbleedmask()
    # app.parse_args(cmdtest.split()[1:])
    # print cmdtest
    # print app(dryrun=True)
    # print

    cmdtest = 'sky_compress --in D00229650_g_03_r0000p01_bleedmasked.fits --skyfilename D00229650_g_03_r0000p01_bleedmask-mini.fits --blocksize 128'
    app = skycompress()
    app.parse_args(cmdtest.split()[1:])
    print cmdtest
    print app(dryrun=True)
    print

    cmdtest = "sky_combine --miniskylist listpcain -o D00229650_g_r0000p01_bleedmask-mini-fp.fits --ccdnums 1,3,4 --invalid S30,N30"

    cmdtest = "sky_fit --infilename D00229650_g_r0000p01_bleedmask-mini-fp.fits --outfilename D00229650_g_r0000p01_skyfit-binned-fp.fits --pcfilename  pca_mini_y2_e2_g_n04.fits"

    cmdtest = "sky_subtract -i  D00229650_g_03_r0000p01_bleedmasked.fits -o D00229650_g_03_r0000p01_skysub.fits --fitfilename D00229650_g_r0000p01_skyfit-binned-fp.fits  --pcfilename skytemp_y2_e2_g_n04_c03.fits --domefilename D_n20150105t0115_g_c03_r2050p02_norm-dflatcor.fits --weight  sky"

    cmdtest = "sex D00229650_g_03_r0000p01_nullwtbkg.fits[0] -c  sexforscamp.config -CHECKIMAGE_TYPE BACKGROUND  -DETECT_THRESH 1000 -FILTER N  -CHECKIMAGE_NAME D00229650_g_03_r0000p01_bkg.fits -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE D00229650_g_03_r0000p01_nullwtbkg.fits[2],D00229650_g_03_r0000p01_nullwtbkg.fits[2] -PARAMETERS_NAME sex.param_bkg -CATALOG_TYPE NONE -INTERP_TYPE ALL -INTERP_MAXXLAG 16 -INTERP_MAXYLAG 16"

    cmdtest = "immask all D00229650_g_03_r0000p01_starflat.fits D00229650_g_03_r0000p01_immask.fits   --minSigma 7.0 --max_angle 75  --max_width 300 --nsig_detect 18 --nsig_mask 12 --nsig_merge 12   --nsig_sky 1.5  --min_fill 0.33  --draw  --write_streaks  --streaksfile D00229650_g_03_r0000p01_streaksfile.fits"

"""
DECam_crosstalk <infile.fits> <outfile> <options>
  -crosstalk <crosstalk matrix- file>
  -linear 
  -crossatthresh <factor> 
  -photflag <0 or 1>
  -presatmask
  -postsatmask
  -overscan
  -overscansample <-1 for MEDIAN, 0 for MEAN, 1 for MEAN w/MINMAX>, 
  -overscanfunction < -50 < N < -1 for cubic SPLINE, 0 for LINE_BY_LINE, 1 < N < 50 for legendre polynomial>
  -overscanorder <1-6>, default=1 <order for legendre polynomial>
  -overscantrim <ncols>, default=0 <trim ncols at both edges of overscan, default = 0>
  -maxhdunum <integer>, default=62 <all valid images>
  -ccdlist <comma-separated list of CCDs to process>
  -hdulist <comma-separated list of HDUs to process>
  -focuschipsout
  -replace <filename>
  -verbose <0-3>
  -help (print help message and exit)
  -version (print version and exit)
"""

"""
usage: pixcorrect_im [-h] [-s SAVECONFIG] [-l LOG] [-v] [-i IN] [-o OUT]
                     [--bias BIAS] [--lincor LINCOR] [--bf BF] [--gain]
                     [--bpm BPM] [--flat FLAT] [--fixcols] [--mini MINI]
                     [--blocksize BLOCKSIZE] [--sky SKY] [--skyfit SKYFIT]
                     [--starflat STARFLAT] [--addweight] [--resaturate]
                     [--null_mask NULL_MASK]
                     [config]

Do image-by-image pixel level corrections

positional arguments:
  config                Configuration file filename

optional arguments:
  -h, --help            show this help message and exit
  -s SAVECONFIG, --saveconfig SAVECONFIG
                        output config file
  -l LOG, --log LOG     the name of the logfile
  -v, --verbose         be verbose
  -i IN, --in IN        input image file name
  -o OUT, --out OUT     output image file name
  --bias BIAS           Bias correction image
  --lincor LINCOR       linearity correction Table
  --bf BF               brighter/fatter correction Table
  --gain                convert ADU to e- using gain values in hdr
  --bpm BPM             bad pixel mask filename
  --flat FLAT           Dome flat correction image
  --fixcols             fix bad columns
  --mini MINI           compressed sky image filename
  --blocksize BLOCKSIZE
                        blocksize for compressed sky image
  --sky SKY             Template file for sky subtraction and weight creation.
                        Requires flat and skyfit files be given.
  --skyfit SKYFIT       MiniDECam file holding sky fit coefficients
  --starflat STARFLAT   Star flat correction image
  --addweight           Add a weight map to the image if none exists
  --resaturate          Put saturated value in BADPIX_SATURATE pixels
  --null_mask NULL_MASK
                        Names of mask bits to null (or an integer mask)
"""


"""
usage: sky_compress [-h] [-s SAVECONFIG] [-l LOG] [-v] [-i IN] [-o OUT]
                    [--skyfilename SKYFILENAME] [--blocksize BLOCKSIZE]
                    [--bitmask BITMASK]
                    [config]

Produce compressed image of sky background

positional arguments:
  config                Configuration file filename

optional arguments:
  -h, --help            show this help message and exit
  -s SAVECONFIG, --saveconfig SAVECONFIG
                        output config file
  -l LOG, --log LOG     the name of the logfile
  -v, --verbose         be verbose
  -i IN, --in IN        input image file name
  -o OUT, --out OUT     output image file name
  --skyfilename SKYFILENAME
                        Filename for compressed sky image
  --blocksize BLOCKSIZE
                        Size of squares in which median is taken for sky
  --bitmask BITMASK     Mask image bits for pixels to ignore in sky estimate
"""

"""
usage: sky_combine [-h] [-s SAVECONFIG] [-l LOG] [-v] [--ccdnums CCDNUMS]
                   [--miniskyfiles MINISKYFILES] [--miniskylist MINISKYLIST]
                   [-o OUTFILENAME] [--mask_value MASK_VALUE]
                   [--invalid INVALID]
                   [config]

Combine sky images of all CCDs in one exposure

positional arguments:
  config                Configuration file filename

optional arguments:
  -h, --help            show this help message and exit
  -s SAVECONFIG, --saveconfig SAVECONFIG
                        output config file
  -l LOG, --log LOG     the name of the logfile
  -v, --verbose         be verbose
  --ccdnums CCDNUMS     Range(s) of ccdnums to combine
  --miniskyfiles MINISKYFILES
                        Filename template for single-chip minisky images
  --miniskylist MINISKYLIST
                        File containing a list of single-chip minisky images
  -o OUTFILENAME, --outfilename OUTFILENAME
                        Filename for combined minisky FITS image
  --mask_value MASK_VALUE
                        Value of pixels without valid sky information
  --invalid INVALID     Value(s) of DETPOS to ignore in sky image

"""

"""
usage: sky_fit [-h] [-s SAVECONFIG] [-l LOG] [-v] [-i IN] [-o OUT]
               [--infilename INFILENAME] [--outfilename OUTFILENAME]
               [--pcfilename PCFILENAME] [--clip_sigma CLIP_SIGMA]
               [config]

Fit coefficients of sky templates to mini-sky image

positional arguments:
  config                Configuration file filename

optional arguments:
  -h, --help            show this help message and exit
  -s SAVECONFIG, --saveconfig SAVECONFIG
                        output config file
  -l LOG, --log LOG     the name of the logfile
  -v, --verbose         be verbose
  -i IN, --in IN        input image file name
  -o OUT, --out OUT     output image file name
  --infilename INFILENAME
                        Filename for input minisky FITS image to fit
  --outfilename OUTFILENAME
                        Filename for minisky FITS image with fit
                        results/resids
  --pcfilename PCFILENAME
                        Filename for minisky principal components
  --clip_sigma CLIP_SIGMA
                        Rejection threshold for robust fitting/statistics

"""

"""
usage: sky_subtract [-h] [-s SAVECONFIG] [-l LOG] [-v] [-i IN] [-o OUT]
                    [--fitfilename FITFILENAME] [--pcfilename PCFILENAME]
                    [--domefilename DOMEFILENAME] [--weight {sky,all,none}]
                    [--resaturate] [--null_mask NULL_MASK]
                    [config]

Subtract sky from images based on principal-component fit and calculate weight
image

positional arguments:
  config                Configuration file filename

optional arguments:
  -h, --help            show this help message and exit
  -s SAVECONFIG, --saveconfig SAVECONFIG
                        output config file
  -l LOG, --log LOG     the name of the logfile
  -v, --verbose         be verbose
  -i IN, --in IN        input image file name
  -o OUT, --out OUT     output image file name
  --fitfilename FITFILENAME
                        Filename for minisky FITS image with PC coefficients
  --pcfilename PCFILENAME
                        Filename for full-res sky principal components
  --domefilename DOMEFILENAME
                        Filename for dome flat (for weight calculation)
  --weight {sky,all,none}
                        Construct weight from sky photons, from all photons,
                        or not at all
  --resaturate          Put saturated value in BADPIX_SATURATE pixels
  --null_mask NULL_MASK
                        Names of mask bits to null (or an integer mask)

"""
"""
SYNTAX: scamp catalog1 [catalog2,...][@catalog_list1 [@catalog_list2 ...]]
[-c <config_file>][-<keyword> <value>]
> to dump a default configuration file: SCAMP -d 
> to dump a default extended configuration file: SCAMP -dd 

PLPLOT-specific options:

Usage:
        scamp [options]

PLplot options:
    -h                   Print out this message
    -v                   Print out the PLplot library version number
    -verbose             Be more verbose than usual
    -debug               Print debugging info (implies -verbose)
    -dev name            Output device name
    -o name              Output filename
    -display name        X server to contact
    -px number           Plots per page in x
    -py number           Plots per page in y
    -geometry geom       Window size/position specified as in X, e.g., 400x300, 400x3
00-100+200, +100-200, etc.
    -wplt xl,yl,xr,yr    Relative coordinates [0-1] of window into plot
    -mar margin          Margin space in relative coordinates (0 to 0.5, def 0)
    -a aspect            Page aspect ratio (def: same as output device)
    -jx justx            Page justification in x (-0.5 to 0.5, def 0)
    -jy justy            Page justification in y (-0.5 to 0.5, def 0)
    -ori orient          Plot orientation (0,1,2,3=landscape,portrait,seascape,upside
-down)
    -freeaspect          Allow aspect ratio to adjust to orientation swaps
    -portrait            Sets portrait mode (both orientation and aspect ratio)
    -width width         Sets pen width (0 <= width)
    -bg color            Background color (FF0000=opaque red, 0000FF_0.1=blue with al
pha of 0.1)
    -ncol0 n             Number of colors to allocate in cmap 0 (upper bound)
    -ncol1 n             Number of colors to allocate in cmap 1 (upper bound)
    -fam                 Create a family of output files
    -fsiz size[kKmMgG]   Output family file size (e.g. -fsiz 0.5G, def MB)
    -fbeg number         First family member number on output
    -finc number         Increment between family members
    -fflen length        Family member number minimum field width
    -nopixmap            Don't use pixmaps in X-based drivers
    -db                  Double buffer X window output
    -np                  No pause between pages
    -server_name name    Main window name of PLplot server (tk driver)
    -dpi dpi             Resolution, in dots per inch (e.g. -dpi 360x360)
    -compression num     Sets compression level in supporting devices
    -cmap0 file name     Initializes color table 0 from a cmap0.pal format file in on
e of standard PLplot paths.
    -cmap1 file name     Initializes color table 1 from a cmap1.pal format file in on
e of standard PLplot paths.
    -locale              Use locale environment (e.g., LC_ALL, LC_NUMERIC, or LANG) t
o set LC_NUMERIC locale (which affects decimal point separator).
    -eofill              For the case where the boundary of the filled region is self
-intersecting, use the even-odd fill rule rather than the default nonzero fill rule.
    -drvopt option[=value][,option[=value]]* Driver specific options

All parameters must be white-space delimited.  Some options are driver
dependent.  Please see the PLplot reference document for more detail.

"""

"""
> SYNTAX: sex <image> [<image2>][-c <configuration_file>][-<keyword> <value>]
> to dump a default configuration file:          sex -d 
> to dump a default extended configuration file: sex -dd 
> to dump a full list of measurement parameters: sex -dp 

"""

"""
    _description = 'Python wrapper around DECam_crosstalk'
    _arguments = (
        (('infile',),
         dict(metavar='infile.fits',help='input file')),
        (('outfile',),
         dict(metavar='outfile.fits',help='output file')),
        (('-crosstalk',),
         dict(metavar='matrix',type=str,help='crosstalk matrix file')),
        (('-linear',),
         dict(action='store_true',help='')),
        (('-crossatthresh',),
          dict(metavar='<factor>',type=float,help='')),
        (('-photflag',),
         dict(type=int,choices=[0,1],help='')),
        (('-presatmask',),
         dict(action='store_true',help='')),
        (('-postsatmask',),
         dict(action='store_true',help='')),
        (('-overscan',),
         dict(action='store_true', help='')),
        (('-overscansample',),
         dict(metavar='type {-1,0,1}', type=int,choices=[-1,0,1],
              help='-1 for MEDIAN, 0 for MEAN, 1 for MEAN w/MINMAX')),
        (('-overscanfunction',),
         dict(metavar='func {-50,50}',default=0,type=int,choices=range(-50,51),
              help='-50 < N < -1 for cubic SPLINE, 0 for LINE_BY_LINE, 1 < N < 50 for legendre polynomial')),
        (('-overscanorder',),
         dict(default=1,type=int,choices=range(1,7),
              help='order for legendre polynomial')),
        (('-overscantrim',),
         dict(metavar='ncols',type=int,default=0,
              help='trim ncols at both edges of overscan')),
        (('-maxhdunum',),
         dict(type=int,default=62, help='all valid images (1-indexed)')),
        (('-ccdlist',),
         dict(help='comma-separated list of CCDs to process')),
        (('-hdulist',),
         dict(help='comma-separated list of HDUs to process')),
        (('-focuschipsout',),
         dict(action='store_true',help='')),
        (('-replace',),
         dict(metavar='filename', help='replace filename')),
        (('-verbose',),
         dict(choices=range(4),default=0,type=int, help='Output verbosity')),
        (('-version',),
         dict(action='store_true', help='print version and exit')),
        (('-help',),
         dict(action='store_true',help='print help message and exit')),
    )
"""

"""
# All of the pixcorrect steps could be done better...
class pixcorrect(DESDMApp):
    _prog = 'pixcorrect_im'
    _section = 'pixcorrect'

    def _create_parser(self):
        from pixcorrect.pixcorrect_im import PixCorrectIm
        parser = PixCorrectIm.parser()        
        parser.prog = self._prog
        parser.description = self._description
        self.namespace.help = '==SUPPRESS=='
        return parser


class skycompress(DESDMApp):
    _prog = 'sky_compress'
    _section = 'skycompress'

    def _create_parser(self):
        from pixcorrect.sky_compress import sky_compress
        parser = sky_compress.parser()        
        parser.prog = self._prog
        parser.description = self._description
        self.namespace.help = '==SUPPRESS=='
        return parser

class skycompress(DESDMApp):
    _prog = 'sky_compress'
    _section = 'skycompress'

    def _create_parser(self):
        from pixcorrect.sky_compress import sky_compress
        parser = sky_compress.parser()        
        parser.prog = self._prog
        parser.description = self._description
        self.namespace.help = '==SUPPRESS=='
        return parser


class skycombine(DESDMApp):
    _prog = 'sky_combine'
    _section = 'skycombine'

    def _create_parser(self):
        from pixcorrect.sky_combine import sky_combine
        parser = sky_compress.parser()        
        parser.prog = self._prog
        parser.description = self._description
        self.namespace.help = '==SUPPRESS=='
        return parser

"""

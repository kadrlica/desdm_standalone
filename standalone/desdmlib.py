#!/usr/bin/env python
"""
Module for wrapping DESDM command line scripts. Most of the work is in
parsing command line arguments into members of the pythonized objects.
"""

import os
import argparse
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import subprocess
from subprocess import Popen,PIPE
from configparser import SafeConfigParser

# Option key string (lowercase)
OPTIONS = '<options>'

class DESDMApp(object):
    """
    Base class for wrapping a DESDM command line script.
    """
    _prog = "desdm_app"
    _description = 'Python wrapper around %s.'
    _section = 'desdmapp'
    _help = '%s -h'
    _cache = None

    def __init__(self,options=None,**kwargs):
        self.section = self._section
        self.namespace = argparse.Namespace()
        self.parser = self._create_parser()
        self._setup_parser()
        if options is None: options = ''
        self.namespace = self.parser.parse_args(options,self.namespace)
        self.namespace.__dict__.update(**kwargs)
        #print '\n',self.__class__.__name__,'\n',self.namespace

    def __call__(self,dryrun=False):
        return self.run(dryrun=dryrun)

    def _create_parser(self):
        prog = self._prog
        description = self._description%self._prog
        parser = ArgumentParser(prog=prog,description=description,
                                formatter_class=ArgumentDefaultsHelpFormatter,
                                add_help=False)
        return parser

    def _setup_parser(self):
        self.parser.prog = self._prog
        self.parser.description = self._description%self._prog

        # The help parser doesn't get added to the namespace
        if self.parser.add_help: 
            self.namespace.help = self.parser.get_default('help')
        # Need to be able to parser without arguments
        for action in self.parser._actions:
            action.required = False

    @classmethod
    def _clear_cache(cls):
        cls._cache = None
        return cls._cache

    def run(self,dryrun=False,verbose=False):
        #self.get_files(verbose=verbose,dryrun=dryrun)
        cmdline = self.cmdline()
        if verbose: 
            name = self.__class__.__name__
            header = "\n{0:-^{width}}\n{1:-^{width}}\n{2:-^{width}}"
            header = header.format(''," Running '%s' "%name,'',width=50)
            print header+'\n'+cmdline+'\n'

        if dryrun: return cmdline

        p = Popen(cmdline,shell=True,stdout=PIPE,stderr=subprocess.STDOUT)
        stdout = []
        for line in iter(p.stdout.readline, b''):
            print line.rstrip()
            stdout.append(line)
        return '\n'.join(stdout)

    def parse_args(self,args=None):
        namespace = getattr(self,'namespace',None)
        if isinstance(args,basestring): args = args.split()
        self.namespace = self.parser.parse_args(args,namespace)
        return self.namespace
        
    def parse_config(self,config,section=None):
        if section is None: section = self.section

        if isinstance(config,basestring):
            configfile = config
            config = SafeConfigParser()
            config.read(configfile)

        # Set the options first (otherwise optional positional arguments disappear)
        if OPTIONS in config._sections[section].keys():
            self.parse_args(config.get(section,OPTIONS))

        # Don't grab the default keys
        dests = [a.dest for a in self.parser._actions]
        for key in config._sections[section].keys():
            if key == '__name__': continue
            if key == OPTIONS: continue
            try: 
                index = dests.index(key)
                action = self.parser._actions[index]
            except ValueError:
                msg = "Unrecognized argument in config: %s\n"%key
                msg += self.get_help()
                raise ValueError(msg)

            if isinstance(action,argparse._StoreTrueAction):
                value = config.getboolean(section,key)
            if isinstance(action,argparse._StoreFalseAction):
                value = config.getboolean(section,key)
            if isinstance(action,argparse._CountAction):
                try:
                    value = int(config.getboolean(section,key))
                except ValueError:
                    value = config.getint(section,key)
            else:
                value = config.get(section,key)
            setattr(self.namespace,key,value)

        return self.namespace

    @classmethod
    def get_help(cls):
        if cls._cache is None:
            cmd = cls._help%cls._prog
            output = subprocess.Popen(cmd.split(),
                                      stdout=PIPE,stderr=subprocess.STDOUT)
            out = output.stdout.read()
            cls._cache = out
        else:
            out = cls._cache
        return out

    @classmethod
    def print_help(cls):
        print cls.get_help()

    def get_files(self,force=False,verbose=False,dryrun=False):
        paths = []
        for key,value in vars(self.namespace).items():
            if not isinstance(value,basestring): continue
            if os.path.exists(value): continue
            if not value.startswith('/pnfs/des/'): continue
            
            path = value
            basename = os.path.basename(path)
            if os.path.exists(basename) and not force: continue
                
            paths += [path]
            if not dryrun: setattr(self.namespace,key,basename)

        if paths:
            import dcache
            dcache.copy(paths,verbose=verbose,dryrun=dryrun)
        return paths

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
    _help = "%s -help"

    def _create_parser(self):
        """
        Create an ArgumentParser from help message.
        """
        parser = ArgumentParser(prog=self._prog,description=self._description,
                                formatter_class=ArgumentDefaultsHelpFormatter,
                                add_help=False)

        try:
            out = self.get_help()
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

class bleedmask(DESDMApp):
    _prog = 'mkbleedmask'
    _section = 'bleedmask'

    def _create_parser(self):
        parser = super(bleedmask,self)._create_parser()
        out = self.get_help()

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

### ASTROMATIC TOOLS ###

class astromatic(DESDMApp):
    _help = '%s -dd'

    def _create_parser(self):
        """
        Create an ArgumentParser.
        """
        parser = super(astromatic,self)._create_parser()

        out = self.get_help()
        lines = map(str.strip,out.strip().split('\n'))
        options = [l for l in lines if (l and not l.startswith('#'))]
        arguments = ['<infile>','<infile2>']

        for arg in arguments:
            args = (arg.strip('<>'),)
            kwargs = dict(nargs='?')
            parser.add_argument(*args,**kwargs)

        # Add the configuration file
        parser.add_argument('-c','--config',help='Configuration file')
        for opt in options:
            olist = opt.lower().split()
            args = ('-'+olist[0],)
            kwargs = dict(action='store',
                          help=opt.split('#')[-1])
            parser.add_argument(*args,**kwargs)

        return parser

class sextractor(astromatic):
    _prog = 'sex'
    _section = 'sextractor'

class psfex(astromatic):
    _prog = 'psfex'
    _section = 'psfex'

class scamp(astromatic):
    _prog = 'scamp'
    _section = 'scamp'

##############################

### Python apps where the parser is inaccessible ###

class PythonApp(DESDMApp):
    
    def _create_parser(self):
        parser = super(PythonApp,self)._create_parser()

        out = self.get_help()
        lines = out.strip().split('\n')
        lines = [l.strip() for l in lines if (len(l)-len(l.lstrip())==2)]

        arguments = [l for l in lines if not l.startswith('-')]
        options = [l for l in lines if l.startswith('-')]

        for arg in arguments:
            args = (arg.split()[0],)
            kwargs = dict(nargs='?')
            parser.add_argument(*args,**kwargs)
         
        # Add the configuration file
        for opt in [o.split(2*' ')[0] for o in options]:
            olist = [o.split() for o in opt.split(',')]
            if len(olist[0]) == 1: action = 'store_true'
            else: action = 'store'
            args = tuple([o[0] for o in olist if o[0].startswith('-')])
            kwargs = dict(action=action)
            parser.add_argument(*args,**kwargs)

        return parser

class fitscombine(PythonApp):
    _prog = 'fitscombine.py'
    _section = 'fitscombine'

class catcombine(PythonApp):
    _prog = 'combine_cats.py'
    _section = 'catcombine'

class updatewcs(PythonApp):
    _prog = 'updateWCS'
    _section = 'updatewcs'

class immask(PythonApp):
    _prog = 'immask all'
    _section = 'immask'


### Pixcorrect functions

class pixcorrect(DESDMApp):
    _prog = 'pixcorrect_im'
    _section = 'pixcorrect_im'

    def _create_parser(self):
        from pixcorrect.pixcorrect_im import PixCorrectIm
        parser = PixCorrectIm.parser()        
        return parser

class nullweight(pixcorrect):
    _prog = 'null_weights'
    _section = 'nullweight'

    def _create_parser(self):
        from pixcorrect.null_weights import null_weights
        parser = null_weights.parser()        
        return parser

class skycompress(DESDMApp):
    _prog = 'sky_compress'
    _section = 'skycompress'
    
    def _create_parser(self):
        from pixcorrect.sky_compress import sky_compress
        parser = sky_compress.parser()        
        return parser

class skycombine(DESDMApp):
    _prog = 'sky_combine'
    _section = 'skycombine'
    
    def _create_parser(self):
        from pixcorrect.sky_combine import sky_combine
        parser = sky_combine.parser()        
        return parser
    
class skyfit(DESDMApp):
    _prog = 'sky_fit'
    _section = 'skyfit'

    def _create_parser(self):
        from pixcorrect.sky_fit import sky_fit
        parser = sky_fit.parser()        
        return parser

class skysubtract(DESDMApp):
    _prog = 'sky_subtract'
    _section = 'skysubtract'

    def _create_parser(self):
        from pixcorrect.sky_subtract import sky_subtract
        parser = sky_subtract.parser()
        return parser
    
class rowinterp(DESDMApp):
    _prog = 'rowinterp_nullweight'
    _section = 'rowinterp'

    def _create_parser(self):
        from pixcorrect.rowinterp_nullweight import RowInterpNullWeight
        parser = RowInterpNullWeight.parser()
        return parser

#--------------------------------#
#      creating regions          #
#--------------------------------#

def read_geometry(infile, outfile, CCD):
    cat = pyfits.open(infile)
    XWIN_IMAGE = cat[2].data.field('XWIN_IMAGE').copy()
    YWIN_IMAGE = cat[2].data.field('YWIN_IMAGE').copy()
    A_IMAGE = cat[2].data.field('A_IMAGE').copy()
    B_IMAGE = cat[2].data.field('B_IMAGE').copy()
    THETA_IMAGE = cat[2].data.field('THETA_IMAGE').copy()
    SPREAD_MODEL = cat[2].data.field('SPREAD_MODEL').copy()
    FLAGS = cat[2].data.field('FLAGS').copy()
    IMAFLAGS_ISO = cat[2].data.field('IMAFLAGS_ISO').copy() 
    del cat
 
    c = 5
    A_IMAGE = c*A_IMAGE
    B_IMAGE = c*B_IMAGE

    color = []
    ellipse = []
    widht = []   


    N = len(SPREAD_MODEL)
    for i in range(N):
        ellipse.append("ellipse")
	if FLAGS[i]<4:
	    if IMAFLAGS_ISO[i] == 0:
	        if SPREAD_MODEL[i] < 0.003 and SPREAD_MODEL[i] > -0.003:	
                    color.append("#color=red") #stars
		    widht.append('width = 4')
		elif SPREAD_MODEL[i] >= 0.003:
		    color.append("#color=yellow") #galaxies
		    widht.append('width = 4')		
	        else:
	            color.append("#color=green") #garbage
		    widht.append('width = 2')
	    elif IMAFLAGS_ISO[i]  & 32 !=0:
		    color.append("#color=blue") #detection near a bright stars.
		    widht.append('width = 3')
	    else:
		color.append("#color=cyan") #object with questionable flags  Flags >=4
                widht.append('width = 3')
	else: 
	    color.append("#color=magenta") #object with questionable flags  Flags >=4
	    widht.append('width = 2')

    DAT =  np.column_stack((ellipse, XWIN_IMAGE, YWIN_IMAGE, A_IMAGE, B_IMAGE, THETA_IMAGE, color, widht))
    np.savetxt(outfile, DAT, delimiter=" ", fmt="%s")   


    
if __name__ == '__main__':

    ### Move this to a tests directory

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


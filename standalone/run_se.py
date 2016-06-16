from desSELib import * 
import numpy as np
import os

#just testing in  fews ccds
EXPNUM =  ConfigSectionMap("General")['expnum']
FILTER = ConfigSectionMap("General")['filter']
NITE = ConfigSectionMap("General")['nite']
CCD = (ConfigSectionMap("General")['chiplist']).split( ',')
rRun = ConfigSectionMap("General")['r']
pRun = ConfigSectionMap("General")['p']
YEAR = ConfigSectionMap("General")['year']
EPOCH = ConfigSectionMap("General")['epoch']

EXPFILE =  'DECam_00'+str(EXPNUM)+'.fits.fz'


args = {'expnum': EXPNUM, 'filter': FILTER, 'ccd':'0', 'r':rRun, 'p':pRun, 'year': YEAR, 'epoch': EPOCH}
#running crosstalk
crosstalk(EXPFILE,NITE,**args)

#running pixelcorrect ans bleedmask
copy_from_Dcash(data_conf+'default.psf')
for ccd in CCD:
        ccdstring= "%02d"%int(ccd)
        pixcorrect(ccdstring,**args)
        sextractor('detrend', ccdstring, **args)
skyCombineFit('mini')

combineFiles('sextractor.fits', 'Scamp_allCCD.fits')
scamp('Scamp_allCCD.fits')
change_head('Scamp_allCCD.head', 'sextractor', 'detrend', CCD, **args)

for ccd in CCD:
        ccdstring= "%02d"%int(ccd)
        bleedmask(ccdstring,'wcs',**args)
        immask(ccdstring,'bleedmasked',**args)

# copy to  Dcash
data_exp = ConfigSectionMap("General")['exp_dir']
dir_nite = data_exp+NITE
os.system('ifdh mkdir '+str(dir_nite) )

dir_final = dir_nite+'/'+EXPNUM
os.system('ifdh mkdir '+str(dir_final) )

os.system('ifdh  cp -D  *sextractor* *immask*  ' +str(dir_final))

import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu


base='/global/homes/z/zonca/m/maps/out'
baseddx9='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/Single_Radiometer'

fwhmdeg=10.
for s in ['1','2','3','4','5']:
    for q in ['18S','23S','19S','22S','20S','21S','18M','23M','19M','22M','20M','21M']:
        fldelta=glob(baseddx9+'/*'+q+'*survey_'+s+'*1s.fits')
        flosgn=glob(base+'/map*osgn*'+q+'*survey'+s+'.fits')
        diffi=pu.get_surv_diff_fname(fldelta[0],flosgn[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diffi,title='SS'+s+' DeltaDX9 - OSGN quad'+q+' smooth 10 deg')
        plt.savefig('plots/osgn_delta/ddx9_osgn_'+q+'_ss'+s+'.png')
        

        
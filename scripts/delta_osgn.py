import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu


base='/global/homes/z/zonca/m/maps/out'


fwhmdeg=10.
for s in ['1','2','3','4','5']:
    for q in ['18','23','19','22','20','21']:
        fldelta=glob(base+'/map*ddx9*'+q+'*survey'+s+'.fits')
        flosgn=glob(base+'/map*osgn*'+q+'*survey'+s+'.fits')
        diffi=pu.get_surv_diff_fname(fldelta[0],flosgn[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diffi,title='SS'+s+' DeltaDX9 - OSGN quad'+Q+' smooth 10 deg')
        plt.savefig('plots/osgn_delta/ddx9_osgn_'+q+'i_ss'+s+'.png')
        
        diffq=pu.get_surv_diff_fname(fldelta[0],flosgn[0],freq=70,fwhm=fwhmdeg,pol='Q')
        hp.mollview(diffi,title='SS'+s+' DeltaDX9 - OSGN quad'+q+' smooth 10 deg')
        plt.savefig('plots/osgn_delta/ddx9_osgn_'+q+'q_ss'+s+'.png')

        diffu=pu.get_surv_diff_fname(fldelta[0],flosgn[0],freq=70,fwhm=fwhmdeg,pol='U')
        hp.mollview(diffi,title='SS'+s+' DeltaDX9 - OSGN quad'+q+' smooth 10 deg')
        plt.savefig('plots/osgn_delta/ddx9_osgn_'+q+'u_ss'+s+'.png')
        
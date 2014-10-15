import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu
basedelta='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/Single_Radiometer'
base='/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/channels'

fwhmdeg=10.
for s in ['1','2','3','4','5']:
    for c in ['18','19','20','21','22','23']:
        fm=glob(basedelta+'/*'+c+'M*'+s+'*.fits')
        fs=glob(basedelta+'/*'+c+'S*'+s+'*.fits')
        diff=pu.get_surv_diff_fname(fm[0],fs[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diff,title='SS'+s+' Delta_DX9 LFI'+c+'M-S  smooth 10 deg')
        plt.savefig('plots/dx9_msdiffs/deltadx9_lfi'+c+'m-s_ss'+s+'.png')
        fm=glob(base+'/*'+c+'M*'+s+'.fits')
        fs=glob(base+'/*'+c+'S*'+s+'.fits')
        diff=pu.get_surv_diff_fname(fm[0],fs[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diff,title='SS'+s+' DX9 LFI'+c+'M-S  smooth 10 deg')
        plt.savefig('plots/dx9_msdiffs/dx9_lfi'+c+'m-s_ss'+s+'.png')

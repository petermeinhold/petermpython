import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu
basedelta='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/Single_Radiometer'
base='/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/channels'
basefull='/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/'
basefulldelta='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
f70delta=glob(basefulldelta+'*70*11*full*.fits')
f70=glob(basefull+'*70*620*full*.fits')

fwhmdeg=10.
for s in ['1','2','3','4','5']:
    #for c in ['18','19','20','21','22','23']:
    for c in ['21','23']:
        fmdelta=glob(basedelta+'/*'+c+'M*'+s+'*.fits')
        fsdelta=glob(basedelta+'/*'+c+'S*'+s+'*.fits')
        diffm=pu.get_surv_diff_fname(fmdelta[0],f70delta[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diffm,title='SS'+s+' Delta_DX9 LFI'+c+'M-freq  smooth 10 deg')
        plt.savefig('plots/dx9_msdiffs/deltadx9_lfi'+c+'m-freq_ss'+s+'.png')
        diffs=pu.get_surv_diff_fname(fsdelta[0],f70delta[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diffs,title='SS'+s+' Delta_DX9 LFI'+c+'S-freq  smooth 10 deg')
        plt.savefig('plots/dx9_msdiffs/deltadx9_lfi'+c+'s-freq_ss'+s+'.png')
        
        fm=glob(base+'/*'+c+'M*'+s+'*.fits')
        fs=glob(base+'/*'+c+'S*'+s+'*.fits')
        diffm=pu.get_surv_diff_fname(fm[0],f70[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diffm,title='SS'+s+' DX9 LFI'+c+'M-freq  smooth 10 deg')
        plt.savefig('plots/dx9_msdiffs/dx9_lfi'+c+'m-freq_ss'+s+'.png')
        diffs=pu.get_surv_diff_fname(fs[0],f70[0],freq=70,fwhm=fwhmdeg)
        hp.mollview(diffs,title='SS'+s+' DX9 LFI'+c+'S-freq  smooth 10 deg')
        plt.savefig('plots/dx9_msdiffs/dx9_lfi'+c+'s-freq_ss'+s+'.png')
        
        
        #fm=glob(base+'/*'+c+'M*'+s+'.fits')
        #fs=glob(base+'/*'+c+'S*'+s+'.fits')
        #diff=pu.get_surv_diff_fname(fm[0],fs[0],freq=70,fwhm=fwhmdeg)
        #hp.mollview(diff,title='SS'+s+' DX9 LFI'+c+'M-S  smooth 10 deg')
        #plt.savefig('plots/dx9_msdiffs/dx9_lfi'+c+'m-s_ss'+s+'.png')

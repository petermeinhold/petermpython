import planck_util as pu
import numpy as np
import cPickle
import pycfitsio
from collections import defaultdict
import healpy as hp
from planck import Planck
from glob import glob
nside=1024
cal='dx9'
chan='30'
base='/global/project/projectdirs/planck/data/mission/DPC_maps/'
dirname=cal
nside=1024
mapname='/lfi/'+'LFI_'+chan+'*'+surv+'.fits'
print base+dirname+mapname
fl30_s1=glob(base+dirname+mapname)
surv='survey_2'
mapname='/lfi/'+'LFI_'+chan+'*'+surv+'.fits'
print base+dirname+mapname
fl30_s2=glob(base+dirname+mapname)
dx9_30_ss1_ss2=pu.get_surv_diff_fname(fl30_s1[0],fl30_s2[0],fwhm=50./60.)


#script to automate masking and rebuilding of full maps tqu
#assume f1 and f2 are full paths to the two maps in question

freq=30
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
m1t=hp.ma(hp.read_map(f1,[0]))
m1q=hp.ma(hp.read_map(f1,[1]))
m1u=hp.ma(hp.read_map(f1,[2]))
m2t=hp.ma(hp.read_map(f2,[0]))
m2q=hp.ma(hp.read_map(f2,[1]))
m2u=hp.ma(hp.read_map(f2,[2]))

tmask=psmask|wmask|(m1t==0).mask|(m2t==0).mask
m1t.mask=tmask
m1q.mask=tmask
m1u.mask=tmask
m2t.mask=tmask
m2q.mask=tmask
m2u.mask=tmask
m1=[m1t.filled(),m1q.filled(),m1u.filled()]
m2=[m2t.filled(),m2q.filled(),m2u.filled()]

cl12tt,cl12ee,cl12bb,cl12te,cl12tb,cl12eb=hp.anafast(m1,m2)
frqac=float(len(m1t[m1t.mask==False])/len(m1t)
l=arange(len(cl12tt)

figure()
labelm1m2='M1 X M2 EE'
plot(l,1e12*l*(l+1)*cl12ee/(frac*2.*np.pi),label=labelm1m2,color='blue')
leg=legend()
pu.thicklegendlines(leg)
title('DX9 30 GHz Ringhalf1, Ringhalf2 cross-spectra (EE) vs Average FSL subtracted version'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


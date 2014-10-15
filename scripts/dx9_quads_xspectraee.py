import planck_util as pu
import numpy as np
import cPickle
import pycfitsio
from collections import defaultdict
import healpy as hp
from planck import Planck
from glob import glob
nside=1024
f1='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/LFI_70_1024_20120611_19_22_nominal.fits'
f2='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/LFI_70_1024_20120611_20_21_nominal.fits'
f3='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/LFI_70_1024_20120617_18_23_nominal.fits'


freq=70
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
m1t=hp.ma(hp.read_map(f1,[0]))
m1q=hp.ma(hp.read_map(f1,[1]))
m1u=hp.ma(hp.read_map(f1,[2]))
m2t=hp.ma(hp.read_map(f2,[0]))
m2q=hp.ma(hp.read_map(f2,[1]))
m2u=hp.ma(hp.read_map(f2,[2]))
m3t=hp.ma(hp.read_map(f3,[0]))
m3q=hp.ma(hp.read_map(f3,[1]))
m3u=hp.ma(hp.read_map(f3,[2]))


tmask=psmask|wmask|(m1t==0).mask|(m2t==0).mask|(m3t==0)
m1t.mask=tmask
m1q.mask=tmask
m1u.mask=tmask
m2t.mask=tmask
m2q.mask=tmask
m2u.mask=tmask
m3t.mask=tmask
m3q.mask=tmask
m3u.mask=tmask
m1=[m1t.filled(),m1q.filled(),m1u.filled()]
m2=[m2t.filled(),m2q.filled(),m2u.filled()]
m3=[m3t.filled(),m3q.filled(),m3u.filled()]

cl12tt,cl12ee,cl12bb,cl12te,cl12tb,cl12eb=hp.anafast(m1,m2)
cl13tt,cl13ee,cl13bb,cl13te,cl13tb,cl13eb=hp.anafast(m1,m3)
cl23tt,cl23ee,cl23bb,cl23te,cl23tb,cl23eb=hp.anafast(m2,m3)
frac=float(len(m1t[m1t.mask==False]))/len(m1t)
l=arange(len(cl12tt))

figure()
labelm1m2='DX9 70 nominal Q1 (1922) x Q2(2021) EE'
plot(l,1e12*l*(l+1)*cl12ee/(frac*2.*np.pi),label=labelm1m2,color='blue')
labelm1m3='DX9 70 nominal Q1 (1922) x Q2(2021) EE'
plot(l,1e12*l*(l+1)*cl13ee/(frac*2.*np.pi),label=labelm1m3,color='green')
labelm2m3='DX9 70 nominal Q1 (2021) x Q2(1823) EE'
plot(l,1e12*l*(l+1)*cl23ee/(frac*2.*np.pi),label=labelm2m3,color='red')
leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 GHz Nominal: Quad cross-spectra (EE)'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


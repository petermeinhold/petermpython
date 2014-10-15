#script to automate masking and rebuilding of full maps tqu
#assume f1 and f2 are full paths to the two maps in question

import planck_util as pu
import numpy as np
import cPickle
import pycfitsio
from collections import defaultdict
import healpy as hp
from planck import Planck
from glob import glob
nside=1024
f1='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/LFI_30_1024_20120611_ringhalf_1_full.fits'
f2='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/LFI_30_1024_20120611_ringhalf_2_full.fits'

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
frac=float(len(m1t[m1t.mask==False]))/len(m1t)
l=arange(len(cl12tt))

figure()
labelm1m2='Ringhalf1 X Ringhalf2 EE'
plot(l,1e12*l*(l+1)*cl12ee/(frac*2.*np.pi),label=labelm1m2,color='blue')
leg=legend()
pu.thicklegendlines(leg)
title('DX9 30 GHz Full: Ringhalf1, Ringhalf2 cross-spectra (EE) vs Average FSL subtracted version'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


#read in the fsl simulations. Assume we can simply add SS1 and SS2 sims in T, Q and U
fslbase='/global/scratch/sd/planck/user/peterm/sidelobes/sims_30ghz_bandpass/average_beam/dx9'
f_fsl1=fslbase+'/madam_slb_30GHz_ss1_512_ave_dx9.fits'
f_fsl2=fslbase+'/madam_slb_30GHz_ss2_512_ave_dx9.fits'

fsl1t=hp.ma(hp.ud_grade(hp.read_map(f_fsl1,[0]),1024))
fsl1q=hp.ma(hp.ud_grade(hp.read_map(f_fsl1,[1]),1024))
fsl1u=hp.ma(hp.ud_grade(hp.read_map(f_fsl1,[2]),1024))

fsl2t=hp.ma(hp.ud_grade(hp.read_map(f_fsl2,[0]),1024))
fsl2q=hp.ma(hp.ud_grade(hp.read_map(f_fsl2,[1]),1024))
fsl2u=hp.ma(hp.ud_grade(hp.read_map(f_fsl2,[2]),1024))

fslat=fsl1t+fsl2t
fslat.mask=fsl1t.mask|fsl2t.mask
fslat.fill_value=0.

fslaq=fsl1q+fsl2q
fslaq.mask=fsl1q.mask|fsl2q.mask
fslaq.fill_value=0.

fslau=fsl1u+fsl2u
fslau.mask=fsl1u.mask|fsl2u.mask
fslau.fill_value=0.

diffm1=[m1t.filled()-fslat.filled(),m1q.filled()-fslaq.filled(),m1u.filled()-fslau.filled()]
diffm2=[m2t.filled()-fslat.filled(),m2q.filled()-fslaq.filled(),m2u.filled()-fslau.filled()]

cl12difftt,cl12diffee,cl12diffbb,cl12diffte,cl12difftb,cl12diffeb=hp.anafast(diffm1,diffm2)
frac=float(len(m1t[m1t.mask==False]))/len(m1t)
l=arange(len(cl12difftt))

labelm1m2='Ringhalf1 - FSL X Ringhalf2 - FSL  EE'
plot(l,1e12*l*(l+1)*cl12diffee/(frac*2.*np.pi),label=labelm1m2,color='red')
leg=legend()
pu.thicklegendlines(leg)
#title('DX9 30 GHz Ringhalf1, Ringhalf2 cross-spectra (EE) vs Average FSL subtracted version'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')
fslat=hp.ma(fslat)
fslaq=hp.ma(fslaq)
fslau=hp.ma(fslau)
fslat.mask=tmask
fslaq.mask=tmask
fslau.mask=tmask
fsla=[fslat.filled(),fslaq.filled(),fslau.filled()]
clfslatt,clfslaee,clfslabb,clfslate,clfslatb,clfslaeb=hp.anafast(fsla)
figure(7)
labelm1m2='FSL sim EE'
plot(l,1e12*l*(l+1)*clfslaee/(frac*2.*np.pi),label=labelm1m2,color='black')
leg=legend()
pu.thicklegendlines(leg)

figure()
labelm1m2='M1 X M2  sim EE'
plot(l,1e12*l*(l+1)*cl12ee/(frac*2.*np.pi),label=labelm1m2,color='blue')
labelm1m2='M1 X M2 -FSL sim EE'
plot(l,1e12*l*(l+1)*cl12diffee/(frac*2.*np.pi),label=labelm1m2,color='red')
labelm1m2='FSL sim EE(auto)'
plot(l,1e12*l*(l+1)*clfslaee/(frac*2.*np.pi),label=labelm1m2,color='magenta')
title('DX9 30 GHz Ringhalf1, Ringhalf2 cross-spectra (EE) vs Average FSL subtracted version'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')
leg=legend()
pu.thicklegendlines(leg)
grid()

figure()
labelm1m2='M1 X M2  sim TT'
plot(l,1e12*l*(l+1)*cl12tt/(frac*2.*np.pi),label=labelm1m2,color='blue')
labelm1m2='M1 X M2 -FSL sim TT'
plot(l,1e12*l*(l+1)*cl12difftt/(frac*2.*np.pi),label=labelm1m2,color='red')
labelm1m2='FSL sim TT(auto)'
plot(l,1e12*l*(l+1)*clfslatt/(frac*2.*np.pi),label=labelm1m2,color='black')


title('DX9 30 GHz Ringhalf1, Ringhalf2 cross-spectra (TT) vs Average FSL subtracted version'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')

leg=legend()
pu.thicklegendlines(leg)
grid()
figure()
labelm1m2='spectrum diff sim EE(auto)'
plot(l,1e12*l*(l+1)*(cl12ee-cl12diffee)/(frac*2.*np.pi),label=labelm1m2,color='blue')

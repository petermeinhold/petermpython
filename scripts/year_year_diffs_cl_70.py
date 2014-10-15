from petermpython import planck_util as pu
import numpy as np
import cPickle
import pycfitsio
from collections import defaultdict
import healpy as hp
from planck import Planck
from glob import glob
nside=1024
freq=70
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
base='/global/project/projectdirs/planck/data/mission/DPC_maps/'

fl=glob(base+'dx9/lfi/LFI_70*7_survey*.fits')
flmr=glob('../gpeterm/maps/DX9_TEST/MultiR/*.fits') 
fl=glob(base+'dx9/lfi/LFI_70*18*_survey*.fits')

l70as1t=hp.ma(hp.read_map(fl[0],0))    
l70as1q=hp.ma(hp.read_map(fl[0],1))    
l70as1u=hp.ma(hp.read_map(fl[0],2))    

l70as2t=hp.ma(hp.read_map(fl[1],0))    
l70as2q=hp.ma(hp.read_map(fl[1],1))    
l70as2u=hp.ma(hp.read_map(fl[1],2))    

l70as3t=hp.ma(hp.read_map(fl[2],0))    
l70as3q=hp.ma(hp.read_map(fl[2],1))    
l70as3u=hp.ma(hp.read_map(fl[2],2))    

l70as4t=hp.ma(hp.read_map(fl[3],0))    
l70as4q=hp.ma(hp.read_map(fl[3],1))    
l70as4u=hp.ma(hp.read_map(fl[3],2))    


tmask=psmask|l70as1t.mask|l70as2t.mask|l70as3t.mask|l70as4t.mask

l70ay1t=(l70as1t+l70as2t)
l70ay1q=(l70as1q+l70as2q)
l70ay1u=(l70as1u+l70as2u)

l70ay2t=(l70as3t+l70as4t)
l70ay2q=(l70as3q+l70as4q)
l70ay2u=(l70as3u+l70as4u)

l70ay1t.mask=tmask
l70ay1q.mask=tmask
l70ay1u.mask=tmask
l70ay1t.mask=tmask
l70ay1q.mask=tmask
l70ay1u.mask=tmask

l70aydiff12=[l70ay1t-l70ay2t,l70ay1q-l70ay2q,l70ay1u-l70ay2u]


frac=float(len(l70ay1q[tmask==False]))/len(tmask)
frac

l70aydiff12allsm=hp.smoothing(l70aydiff12,fwhm=2.*np.pi/180.)

psyyt=pu.rescan_to_psmap(l70aydiff12allsm[0],startring=0,stopring=1100) 
psyyq=pu.rescan_to_psmap(l70aydiff12allsm[1],startring=0,stopring=1100) 
psyyu=pu.rescan_to_psmap(l70aydiff12allsm[2],startring=0,stopring=1100) 


cl70ydifftt,cl70ydiffee,cl70ydiffbb,cl70ydiffte,cl70ydifftb,cl70ydiffeb=hp.anafast(l70aydiff12)
l=arange(len(cl70ydifftt))

figure()
plot l,l*(l+1)*cl70ydifftt/frac,label='DX9 70 Y 1-2 TT'
plot l,l*(l+1)*cl70ydiffee/frac,label='DX9 70 Y 1-2 EE'
plot l,l*(l+1)*cl70ydiffbb/frac,label='DX9 70 Y 1-2 BB'
leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 Y 1-2  Quad 18-23'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


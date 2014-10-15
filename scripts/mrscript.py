import planck_util as pu
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

dx9mr_s1=hp.read_map(flmr[1],[0,1,2]) 
dx9mr_s2=hp.read_map(flmr[0],[0,1,2]) 
dx9mr_s3=hp.read_map(flmr[2],[0,1,2]) 
dx9mr_s4=hp.read_map(flmr[3],[0,1,2]) 

dx9mr_s1t=hp.ma(dx9mr_s1[0])
dx9mr_s1q=hp.ma(dx9mr_s1[1])
dx9mr_s1u=hp.ma(dx9mr_s1[2])

dx9mr_s2t=hp.ma(dx9mr_s2[0])
dx9mr_s2q=hp.ma(dx9mr_s2[1])
dx9mr_s2u=hp.ma(dx9mr_s2[2])

dx9mr_s3t=hp.ma(dx9mr_s3[0])
dx9mr_s3q=hp.ma(dx9mr_s3[1])
dx9mr_s3u=hp.ma(dx9mr_s3[2])

dx9mr_s4t=hp.ma(dx9mr_s4[0])
dx9mr_s4q=hp.ma(dx9mr_s4[1])
dx9mr_s4u=hp.ma(dx9mr_s4[2])

tmaskmr=psmask|wmask|(dx9mr_s1t==0).mask|(dx9mr_s2t==0).mask|(dx9mr_s3t==0).mask|(dx9mr_s4t==0).mask

dx9mr_s1t.mask=tmaskmr
dx9mr_s1q.mask=tmaskmr
dx9mr_s1u.mask=tmaskmr

dx9mr_s2t.mask=tmaskmr
dx9mr_s2q.mask=tmaskmr
dx9mr_s2u.mask=tmaskmr

dx9mr_s3t.mask=tmaskmr
dx9mr_s3q.mask=tmaskmr
dx9mr_s3u.mask=tmaskmr

dx9mr_s4t.mask=tmaskmr
dx9mr_s4q.mask=tmaskmr
dx9mr_s4u.mask=tmaskmr

dx9mr_s1=[dx9mr_s1t.filled(),dx9mr_s1q.filled(),dx9mr_s1u.filled()]
dx9mr_s2=[dx9mr_s2t.filled(),dx9mr_s2q.filled(),dx9mr_s2u.filled()]
dx9mr_s3=[dx9mr_s3t.filled(),dx9mr_s3q.filled(),dx9mr_s3u.filled()]
dx9mr_s4=[dx9mr_s4t.filled(),dx9mr_s4q.filled(),dx9mr_s4u.filled()]

dx9_s1=hp.read_map(fl[0],[0,1,2]) 
dx9_s2=hp.read_map(fl[1],[0,1,2]) 
dx9_s3=hp.read_map(fl[2],[0,1,2]) 
dx9_s4=hp.read_map(fl[3],[0,1,2]) 

dx9_s1t=hp.ma(dx9_s1[0])
dx9_s1q=hp.ma(dx9_s1[1])
dx9_s1u=hp.ma(dx9_s1[2])

dx9_s2t=hp.ma(dx9_s2[0])
dx9_s2q=hp.ma(dx9_s2[1])
dx9_s2u=hp.ma(dx9_s2[2])

dx9_s3t=hp.ma(dx9_s3[0])
dx9_s3q=hp.ma(dx9_s3[1])
dx9_s3u=hp.ma(dx9_s3[2])

dx9_s4t=hp.ma(dx9_s4[0])
dx9_s4q=hp.ma(dx9_s4[1])
dx9_s4u=hp.ma(dx9_s4[2])

tmask=psmask|wmask|(dx9_s1t==0).mask|(dx9_s2t==0).mask|(dx9_s3t==0).mask|(dx9_s4t==0).mask

dx9_s1t.mask=tmask
dx9_s1q.mask=tmask
dx9_s1u.mask=tmask

dx9_s2t.mask=tmask
dx9_s2q.mask=tmask
dx9_s2u.mask=tmask

dx9_s3t.mask=tmask
dx9_s3q.mask=tmask
dx9_s3u.mask=tmask

dx9_s4t.mask=tmask
dx9_s4q.mask=tmask
dx9_s4u.mask=tmask

dx9_s1=[dx9_s1t.filled(),dx9_s1q.filled(),dx9_s1u.filled()]
dx9_s2=[dx9_s2t.filled(),dx9_s2q.filled(),dx9_s2u.filled()]
dx9_s3=[dx9_s3t.filled(),dx9_s3q.filled(),dx9_s3u.filled()]
dx9_s4=[dx9_s4t.filled(),dx9_s4q.filled(),dx9_s4u.filled()]

cl12tt,cl12ee,cl12bb,cl12te,cl12tb,cl12eb=hp.anafast(dx9_s1,dx9_s2)
cl13tt,cl13ee,cl13bb,cl13te,cl13tb,cl13eb=hp.anafast(dx9_s1,dx9_s3)
cl24tt,cl24ee,cl24bb,cl24te,cl24tb,cl24eb=hp.anafast(dx9_s2,dx9_s4)

clmr12tt,clmr12ee,clmr12bb,clmr12te,clmr12tb,clmr12eb=hp.anafast(dx9mr_s1,dx9mr_s2)
clmr13tt,clmr13ee,clmr13bb,clmr13te,clmr13tb,clmr13eb=hp.anafast(dx9mr_s1,dx9mr_s3)
clmr24tt,clmr24ee,clmr24bb,clmr24te,clmr24tb,clmr24eb=hp.anafast(dx9mr_s2,dx9mr_s4)

frac=float(len(dx9_s1t[dx9_s1t.mask == False]))/len(dx9_s1t)
fracmr=float(len(dx9mr_s1t[dx9mr_s1t.mask == False]))/len(dx9mr_s1t)

l=arange(len(clmr12tt))
figure()
plot(l,1e12*l*(l+1)*cl12ee/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS2 EE',color='blue')
plot(l,1e12*l*(l+1)*cl13ee/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS3 EE',color='green')
plot(l,1e12*l*(l+1)*cl24ee/(frac*2.*np.pi),label='DX9 70 GHz SS2 x SS4 EE',color='red')

plot(l,1e12*l*(l+1)*clmr12ee/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS2 EE',color='blue',lw=2)
plot(l,1e12*l*(l+1)*clmr13ee/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS3 EE',color='green',lw=2)
plot(l,1e12*l*(l+1)*clmr24ee/(frac*2.*np.pi),label='DX9 MR 70 GHz SS2 x SS4 EE',color='red',lw=2)
leg=legend()

title('DX9 70 GHz 18+23 survey cross-spectra (EE) vs DX9-multiR'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')

figure()
plot(l,1e12*l*(l+1)*cl12te/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS2 TE',color='blue')
plot(l,1e12*l*(l+1)*cl13te/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS3 TE',color='green')
plot(l,1e12*l*(l+1)*cl24te/(frac*2.*np.pi),label='DX9 70 GHz SS2 x SS4 TE',color='red')

plot(l,1e12*l*(l+1)*clmr12te/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS2 TE',color='blue',lw=2)
plot(l,1e12*l*(l+1)*clmr13te/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS3 TE',color='green',lw=2)
plot(l,1e12*l*(l+1)*clmr24te/(frac*2.*np.pi),label='DX9 MR 70 GHz SS2 x SS4 TE',color='red',lw=2)
leg=legend()

title('DX9 70 GHz survey cross-spectra (TE) vs DX9-multiR'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


figure()
plot(l,1e12*l*(l+1)*cl12tt/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS2 TT',color='blue')
plot(l,1e12*l*(l+1)*cl13tt/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS3 TT',color='green')
plot(l,1e12*l*(l+1)*cl24tt/(frac*2.*np.pi),label='DX9 70 GHz SS2 x SS4 TT',color='red')

plot(l,1e12*l*(l+1)*clmr12tt/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS2 TT',color='blue',lw=2)
plot(l,1e12*l*(l+1)*clmr13tt/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS3 TT',color='green',lw=2)
plot(l,1e12*l*(l+1)*clmr24tt/(frac*2.*np.pi),label='DX9 MR 70 GHz SS2 x SS4 TT',color='red',lw=2)
leg=legend()

title('DX9 70 GHz (18+23)survey cross-spectra (TT) vs DX9-multiR'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')



figure()
plot(l,1e12*l*(l+1)*cl12bb/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS2 BB',color='blue')
plot(l,1e12*l*(l+1)*cl13bb/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS3 BB',color='green')
plot(l,1e12*l*(l+1)*cl24bb/(frac*2.*np.pi),label='DX9 70 GHz SS2 x SS4 BB',color='red')

plot(l,1e12*l*(l+1)*clmr12bb/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS2 BB',color='blue',lw=2)
plot(l,1e12*l*(l+1)*clmr13bb/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS3 BB',color='green',lw=2)
plot(l,1e12*l*(l+1)*clmr24bb/(frac*2.*np.pi),label='DX9 MR 70 GHz SS2 x SS4 BB',color='red',lw=2)
leg=legend()

title('DX9 70 GHz (18+23) survey cross-spectra (BB) vs DX9-multiR'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


figure()
plot(l,1e12*l*(l+1)*cl12bb/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS2 BB',color='blue')
plot(l,1e12*l*(l+1)*cl13bb/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS3 BB',color='blue')
plot(l,1e12*l*(l+1)*cl24bb/(frac*2.*np.pi),label='DX9 70 GHz SS2 x SS4 BB',color='blue')

plot(l,1e12*l*(l+1)*clmr12bb/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS2 BB',color='red')
plot(l,1e12*l*(l+1)*clmr13bb/(frac*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS3 BB',color='red')
plot(l,1e12*l*(l+1)*clmr24bb/(frac*2.*np.pi),label='DX9 MR 70 GHz SS2 x SS4 BB',color='red')
leg=legend()

title('DX9 70 GHz survey cross-spectra (BB) vs DX9-multiR'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')

figure()


plot(l,1e12*l*(l+1)*clmr12ee/(fracmr*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS2 EE',color='red')
plot(l,1e12*l*(l+1)*clmr13ee/(fracmr*2.*np.pi),label='DX9 MR 70 GHz SS1 x SS3 EE',color='red')
plot(l,1e12*l*(l+1)*clmr24ee/(fracmr*2.*np.pi),label='DX9 MR 70 GHz SS2 x SS4 EE',color='red')

plot(l,1e12*l*(l+1)*cl12ee/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS2 EE',color='blue')
plot(l,1e12*l*(l+1)*cl13ee/(frac*2.*np.pi),label='DX9 70 GHz SS1 x SS3 EE',color='blue')
plot(l,1e12*l*(l+1)*cl24ee/(frac*2.*np.pi),label='DX9 70 GHz SS2 x SS4 EE',color='blue')
leg=legend()

title('DX9 70 GHz survey cross-spectra (EE) vs DX9-multiR >>> 18+23 only!'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')





nside=1024
freq=70
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
base='/global/project/projectdirs/planck/data/mission/DPC_maps/'


sl70g_s1t=hp.ma(sl70_s1[0])
sl70g_s1q=hp.ma(sl70_s1[1])
sl70g_s1u=hp.ma(sl70_s1[2])

sl70g_s2t=hp.ma(sl70_s2[0])
sl70g_s2q=hp.ma(sl70_s2[1])
sl70g_s2u=hp.ma(sl70_s2[2])
tmask70=psmask|wmask|(sl70_s1q==0).mask|(sl70_s2q==0).mask

sl70g_s1t.mask=tmask70
sl70g_s1q.mask=tmask70
sl70g_s1u.mask=tmask70

sl70g_s2t.mask=tmask70
sl70g_s2q.mask=tmask70
sl70g_s2u.mask=tmask70

sl70g_s1all=[sl70g_s1t.filled(),sl70g_s1q.filled(),sl70g_s1u.filled()]
sl70g_s2all=[sl70g_s2t.filled(),sl70g_s2q.filled(),sl70g_s2u.filled()]

clsl70gtt,clsl70gee,clsl70gbb,clsl70gte,clsl70gtb,clsl70geb=hp.anafast(sl70g_s1all,sl70g_s2all)

fracsl70g=float(len(sl70g_s1t[sl70g_s1t.mask == False]))/len(sl70g_s1t)
plot(l,1e12*l*(l+1)*clsl70gee/(fracsl70g*2.*np.pi),label='DX9 70 GHz FSL galaxy sim SS1 x SS2 EE')

leg=legend()

title('70 GHz FSL Galaxy sim, survey cross-spectra (EE)'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')



sl70gb_s1t=hp.ma(sl70_s1[0])
sl70gb_s1q=hp.ma(sl70_s1[1])
sl70gb_s1u=hp.ma(sl70_s1[2])

sl70gb_s2t=hp.ma(sl70_s2[0])
sl70gb_s2q=hp.ma(sl70_s2[1])
sl70gb_s2u=hp.ma(sl70_s2[2])
tmask70b=psmask|wmask|(sl70gb_s1q==0).mask|(sl70gb_s2q==0).mask

sl70gb_s1t.mask=tmask70b
sl70gb_s1q.mask=tmask70b
sl70gb_s1u.mask=tmask70b

sl70gb_s2t.mask=tmask70b
sl70gb_s2q.mask=tmask70b
sl70gb_s2u.mask=tmask70b

sl70gb_s1all=[sl70gb_s1t.filled(),sl70gb_s1q.filled(),sl70gb_s1u.filled()]
sl70gb_s2all=[sl70gb_s2t.filled(),sl70gb_s2q.filled(),sl70gb_s2u.filled()]

clsl70gbtt,clsl70gbee,clsl70gbbb,clsl70gbte,clsl70gbtb,clsl70gbeb=hp.anafast(sl70gb_s1all,sl70gb_s2all)
figure()
fracsl70gb=float(len(sl70gb_s1t[sl70gb_s1t.mask == False]))/len(sl70gb_s1t)
plot(l,1e12*l*(l+1)*clsl70gbee/(fracsl70gb*2.*np.pi),label='DX9 70 GHz FSL galaxy sim SS1 x SS2 EE binned')

leg=legend()

title('70 GHz FSL Galaxy sim, survey cross-spectra (EE)'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


#Andrea's tests using HFI destriping masks
flhfi=glob(base+'map*hfidx7*2.fits')
fllfi=glob(base+'map*lfidx9*2.fits')

flhfi
['/global/project/projectdirs/planck/user/zonca/issues/dx9_toast_wikipost/maps/map_dx9_070_nominalhfidx7mask_nofilt_subchunk_1_of_2.fits',
 '/global/project/projectdirs/planck/user/zonca/issues/dx9_toast_wikipost/maps/map_dx9_070_nominalhfidx7mask_nofilt_subchunk_2_of_2.fits']
fllfi
['/global/project/projectdirs/planck/user/zonca/issues/dx9_toast_wikipost/maps/map_dx9_070_nominallfidx9mask_nofilt_subchunk_1_of_2.fits',
 '/global/project/projectdirs/planck/user/zonca/issues/dx9_toast_wikipost/maps/map_dx9_070_nominallfidx9mask_nofilt_subchunk_2_of_2.fits']



freq=70
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))

hfimaskt_1=hp.ma(hp.read_map(flhfi[0],[0]))
hfimaskq_1=hp.ma(hp.read_map(flhfi[0],[1]))
hfimasku_1=hp.ma(hp.read_map(flhfi[0],[2]))

hfimaskt_2=hp.ma(hp.read_map(flhfi[1],[0]))
hfimaskq_2=hp.ma(hp.read_map(flhfi[1],[1]))
hfimasku_2=hp.ma(hp.read_map(flhfi[1],[2]))


dx9maskt_1=hp.ma(hp.read_map(fllfi[0],[0]))
dx9maskq_1=hp.ma(hp.read_map(fllfi[0],[1]))
dx9masku_1=hp.ma(hp.read_map(fllfi[0],[2]))

dx9maskt_2=hp.ma(hp.read_map(fllfi[1],[0]))
dx9maskq_2=hp.ma(hp.read_map(fllfi[1],[1]))
dx9masku_2=hp.ma(hp.read_map(fllfi[1],[2]))



tmask=psmask|wmask|(hfimaskt_1==0).mask|(dx9maskt_1==0).mask|(hfimaskt_2==0).mask|(dx9maskt_2==0).mask

hfimaskt_1.mask=tmask
hfimaskt_1.mask=tmask
hfimaskt_1.mask=tmask
hfimaskt_2.mask=tmask
hfimaskt_2.mask=tmask
hfimaskt_2.mask=tmask

dx9maskt_1.mask=tmask
dx9maskt_1.mask=tmask
dx9maskt_1.mask=tmask
dx9maskt_2.mask=tmask
dx9maskt_2.mask=tmask
dx9maskt_2.mask=tmask

hfimaskall1=[hfimaskt_1.filled(),hfimaskq_1.filled(),hfimasku_1.filled()]
dx9maskall1=[dx9maskt_1.filled(),dx9maskq_1.filled(),dx9masku_1.filled()]

hfimaskall2=[hfimaskt_2.filled(),hfimaskq_2.filled(),hfimasku_2.filled()]
dx9maskall2=[dx9maskt_2.filled(),dx9maskq_2.filled(),dx9masku_2.filled()]

clhfimask12tt,clhfimask12ee,clhfimask12bb,clhfimask12te,clhfimask12tb,clhfimask12eb=hp.anafast(hfimaskall1,hfimaskall2)

cldx9mask12tt,cldx9mask12ee,cldx9mask12bb,cldx9mask12te,cldx9mask12tb,cldx9mask12eb=hp.anafast(dx9maskall1,dx9maskall2)

l=arange(len(clhfimask12tt))

figure()
fracs=float(len(hfimaskt_1[hfimaskt_1.mask == False]))/len(hfimaskt_1)

plot(l,1e12*l*(l+1)*clhfimask12ee/(fracs*2.*np.pi),label='DX9 70 GHz nominal ringhalf 1 x ringhalf 2 EE')
plot(l,1e12*l*(l+1)*cldx9mask12ee/(fracs*2.*np.pi),label='DX9 HFImask 70 GHz ringhalf 1 x ringhalf 2 EE')

leg=legend()

title('70 GHz DX9 nominal vs HFI mask for destriping, first/second halfring cross-spectra (EE)'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')





from petermpython import planck_util as pu

from glob import glob
import healpy as hp

freq=70
nside=1024
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))

dx9base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/'

fl1s=glob(dx9base+'*1s.fits')


fl1m=glob(dx9base+'*70*61*nominal.fits')

m70_1s_at=hp.ma(hp.read_map(fl1s[2],0))
m70_1s_aq=hp.ma(hp.read_map(fl1s[2],1))
m70_1s_au=hp.ma(hp.read_map(fl1s[2],2))

m70_1s_bt=hp.ma(hp.read_map(fl1s[3],0))
m70_1s_bq=hp.ma(hp.read_map(fl1s[3],1))
m70_1s_bu=hp.ma(hp.read_map(fl1s[3],2))

m70_1s_ct=hp.ma(hp.read_map(fl1s[4],0))
m70_1s_cq=hp.ma(hp.read_map(fl1s[4],1))
m70_1s_cu=hp.ma(hp.read_map(fl1s[4],2))

tmask=psmask|wmask|m70_1s_aq.mask|m70_1s_bq.mask|m70_1s_cq.mask
frac1s=len(tmask==False)/len(tmask)
m70_1s_at.mask=tmask
m70_1s_aq.mask=tmask
m70_1s_au.mask=tmask
m70_1s_bt.mask=tmask
m70_1s_bq.mask=tmask
m70_1s_bu.mask=tmask
m70_1s_ct.mask=tmask
m70_1s_cq.mask=tmask
m70_1s_cu.mask=tmask

m70_1s_a=[m70_1s_at,m70_1s_aq,m70_1s_au]
m70_1s_b=[m70_1s_bt,m70_1s_bq,m70_1s_bu]
m70_1s_c=[m70_1s_ct,m70_1s_cq,m70_1s_cu]

clxabtt,clxabee,clxabbb,clxabte,clxabtb,clxabeb=hp.anafast(m70_1s_a,m70_1s_b)
clxactt,clxacee,clxacbb,clxacte,clxactb,clxaceb=hp.anafast(m70_1s_a,m70_1s_c)
clxbctt,clxbcee,clxbcbb,clxbcte,clxbctb,clxbceb=hp.anafast(m70_1s_b,m70_1s_c)

m70_1m_at=hp.ma(hp.read_map(fl1m[2],0))
m70_1m_aq=hp.ma(hp.read_map(fl1m[2],1))
m70_1m_au=hp.ma(hp.read_map(fl1m[2],2))

m70_1m_bt=hp.ma(hp.read_map(fl1m[0],0))
m70_1m_bq=hp.ma(hp.read_map(fl1m[0],1))
m70_1m_bu=hp.ma(hp.read_map(fl1m[0],2))

m70_1m_ct=hp.ma(hp.read_map(fl1m[1],0))
m70_1m_cq=hp.ma(hp.read_map(fl1m[1],1))
m70_1m_cu=hp.ma(hp.read_map(fl1m[1],2))

tmask=psmask|wmask|m70_1m_aq.mask|m70_1m_bq.mask|m70_1m_cq.mask
frac1m=len(tmask==False)/len(tmask)
m70_1m_at.mask=tmask
m70_1m_aq.mask=tmask
m70_1m_au.mask=tmask
m70_1m_bt.mask=tmask
m70_1m_bq.mask=tmask
m70_1m_bu.mask=tmask
m70_1m_ct.mask=tmask
m70_1m_cq.mask=tmask
m70_1m_cu.mask=tmask

m70_1m_a=[m70_1m_at,m70_1m_aq,m70_1m_au]
m70_1m_b=[m70_1m_bt,m70_1m_bq,m70_1m_bu]
m70_1m_c=[m70_1m_ct,m70_1m_cq,m70_1m_cu]

cl1mabtt,cl1mabee,cl1mabbb,cl1mabte,cl1mabtb,cl1mabeb=hp.anafast(m70_1m_a,m70_1m_b)
cl1mactt,cl1macee,cl1macbb,cl1macte,cl1mactb,cl1maceb=hp.anafast(m70_1m_a,m70_1m_c)
cl1mbctt,cl1mbcee,cl1mbcbb,cl1mbcte,cl1mbctb,cl1mbceb=hp.anafast(m70_1m_b,m70_1m_c)
l=arange(len(clxabtt))
figure()
plot l,l*(l+1)*cl1mabee/frac1m,label='70 GHz 1m quads AxB EE'
plot l,l*(l+1)*clxabee/frac1s,label='70 GHz 1s quads AxB EE'

leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 GHz Xspectra 1 minute vs 1 second baselines'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')      

figure()
plot l,l*(l+1)*cl1macee/frac1m,label='70 GHz 1m quads AxC EE'
plot l,l*(l+1)*clxacee/frac1s,label='70 GHz 1s quads AxC EE'
leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 GHz Xspectra 1 minute vs 1 second baselines'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')      

figure()
plot l,l*(l+1)*cl1mbcee/frac1m,label='70 GHz 1m quads BxC EE'
plot l,l*(l+1)*clxbcee/frac1s,label='70 GHz 1s quads BxC EE'
leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 GHz Xspectra 1 minute vs 1 second baselines'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')      



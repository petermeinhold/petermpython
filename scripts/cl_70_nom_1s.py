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

fl=glob(base+'dx9/lfi/LFI_70*nominal.fits')

fl
m1922nom=hp.read_map(fl[0],(0,1,2))

m2021nom=hp.read_map(fl[1],(0,1,2))

m1823nom=hp.read_map(fl[2],(0,1,2))



m1922nomt=hp.ma(m1922nom[0])

m1922nomq=hp.ma(m1922nom[1])

m1922nomu=hp.ma(m1922nom[2])

m1922nomt.mask=psmask|wmask

m1922nomq.mask=psmask|wmask

m1922nomu.mask=psmask|wmask



m1823nomt=hp.ma(m1823nom[0])

m1823nomq=hp.ma(m1823nom[1])

m1823nomu=hp.ma(m1823nom[2])

m1823nomt.mask=psmask|wmask

m1823nomq.mask=psmask|wmask

m1823nomu.mask=psmask|wmask



m2021nomt=hp.ma(m2021nom[0])

m2021nomq=hp.ma(m2021nom[1])

m2021nomu=hp.ma(m2021nom[2])

m2021nomt.mask=psmask|wmask

m2021nomq.mask=psmask|wmask

m2021nomu.mask=psmask|wmask




m1922=[m1922nomt.filled(),m1922nomq.filled(),m1922nomu.filled()]
m1823=[m1823nomt.filled(),m1823nomq.filled(),m1823nomu.filled()]
m2021=[m2021nomt.filled(),m2021nomq.filled(),m2021nomu.filled()]

clabtt,clabee,clabbb,clabte,clabtb,clabeb=hp.anafast(m1922,m1823)

clactt,clacee,clacbb,clacte,clactb,claceb=hp.anafast(m1922,m2021)

clbctt,clbcee,clbcbb,clbcte,clbctb,clbceb=hp.anafast(m1823,m2021)


fl1s=glob(base+'dx9/lfi/LFI_70*1s.fits')

fl1s
m1922nom1s=hp.read_map(fl1s[0],(0,1,2))

m2021nom1s=hp.read_map(fl1s[1],(0,1,2))

m1823nom1s=hp.read_map(fl1s[2],(0,1,2))



m1922nom1st=hp.ma(m1922nom1s[0])

m1922nom1sq=hp.ma(m1922nom1s[1])

m1922nom1su=hp.ma(m1922nom1s[2])

m1922nom1st.mask=psmask|wmask

m1922nom1sq.mask=psmask|wmask

m1922nom1su.mask=psmask|wmask



m1823nom1st=hp.ma(m1823nom1s[0])

m1823nom1sq=hp.ma(m1823nom1s[1])

m1823nom1su=hp.ma(m1823nom1s[2])

m1823nom1st.mask=psmask|wmask

m1823nom1sq.mask=psmask|wmask

m1823nom1su.mask=psmask|wmask



m2021nom1st=hp.ma(m2021nom1s[0])

m2021nom1sq=hp.ma(m2021nom1s[1])

m2021nom1su=hp.ma(m2021nom1s[2])

m2021nom1st.mask=psmask|wmask

m2021nom1sq.mask=psmask|wmask

m2021nom1su.mask=psmask|wmask


m1922_1s=[m1922nom1st.filled(),m1922nom1sq.filled(),m1922nom1su.filled()]
m1823_1s=[m1823nom1st.filled(),m1823nom1sq.filled(),m1823nom1su.filled()]
m2021_1s=[m2021nom1st.filled(),m2021nom1sq.filled(),m2021nom1su.filled()]

cl1sabtt,cl1sabee,cl1sabbb,cl1sabte,cl1sabtb,cl1sabeb=hp.anafast(m1922_1s,m1823_1s)

cl1sactt,cl1sacee,cl1sacbb,cl1sacte,cl1sactb,cl1saceb=hp.anafast(m1922_1s,m2021_1s)

cl1sbctt,cl1sbcee,cl1sbcbb,cl1sbcte,cl1sbctb,cl1sbceb=hp.anafast(m1823_1s,m2021_1s)


frac=float(len(m1922nomt[m1922nomt.mask == False]))/len(m1922nomt)
#estimate sky fraction to correct pseudospectra to first order
plot l,1e12*l*(l+1)*clabee/(frac*2.*np.pi),label='LFI1922 x 1823 60s base'
plot l,1e12*l*(l+1)*cl1sabee/(frac*2.*np.pi),label='LFI1922 x 1823 1s base'
plot l,1e12*l*(l+1)*clacee/(frac*2.*np.pi),label='LFI1922 x 2021 60s base'
plot l,1e12*l*(l+1)*cl1sacee/(frac*2.*np.pi),label='LFI1922 x 2021 1s base'
leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 GHz x spectra (EE) showing low ell excess'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')


dst12i=pu.get_dst_map(surv='surv1-2',cal='DX8S_70',pol='I')
dst12q=pu.get_dst_map(surv='surv1-2',cal='DX8S_70',pol='Q')
dst12u=pu.get_dst_map(surv='surv1-2',cal='DX8S_70',pol='U')
dst34i=pu.get_dst_map(surv='surv3-4',cal='DX8S_70',pol='I')
dst34q=pu.get_dst_map(surv='surv3-4',cal='DX8S_70',pol='Q')
dst34u=pu.get_dst_map(surv='surv3-4',cal='DX8S_70',pol='U')


dst12i.mask=psmask|wmask
dst12q.mask=psmask|wmask
dst12u.mask=psmask|wmask
dst34i.mask=psmask|wmask
dst34q.mask=psmask|wmask
dst34u.mask=psmask|wmask

dst12=[dst12i.filled(),dst12q.filled(),dst12u.filled()]
dst34=[dst34i.filled(),dst34q.filled(),dst34u.filled()]
cldstabtt,cldstabee,cldstabbb,cldstabte,cldstabtb,cldstabeb=hp.anafast(dst12,dst34)



fl70=glob('../fsl_70ghz_sim/*.fits')
fl70
#  ['../fsl_70ghz_sim/madam_fgslb_70GHz_all_ss2_1024binmap.fits', '../fsl_70ghz_sim/madam_fgslb_70GHz_all_ss1_1024binmap.fits']

fsl70ss1=hp.read_map(fl70[0],(0,1,2))
fsl70ss2=hp.read_map(fl70[0],(0,1,2))
fsl70ss1=hp.read_map(fl70[1],(0,1,2))
fsl70i=hp.ma(fsl70ss1[0]+fsl70ss2[0])/2.
fsl70q=hp.ma(fsl70ss1[1]+fsl70ss2[1])/2.
fsl70u=hp.ma(fsl70ss1[2]+fsl70ss2[2])/2.
fsl70i.mask=psmask|wmask
fsl70q.mask=psmask|wmask
fsl70u.mask=psmask|wmask

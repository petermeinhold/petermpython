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

nside=128
nphi=500
stepring=10
vecs=qr.get_hfi_ring_spin_axes()* np.cos( ang * np.pi / 180. )

startring=0
stopring=5500
x=arange(nphi)
psmap=[]
for ring in range(startring,stopring,stepring):
    vec = vecs[ring] 
    xs, ys, zs = qr.maths.get_ring_for_spin_vec(vec, npts=nphi)
    tht   = np.pi*0.5 - np.arcsin(zs) #latitude
    phi   = np.arctan2(ys, xs)        #longitude
    thtpneg=tht[phi<0]
    thtmin=x[(tht-np.pi/2)==np.min(thtpneg-np.pi/2)]
    tpix=np.roll(hp.ang2pix(nside,tht,phi),thtmin[0])
    psmap.append(dx930[tpix])

startring=5500
stopring=10500
psmap2=[]
for ring in range(startring,stopring,stepring):
    vec = vecs[ring] 
    xs, ys, zs = qr.maths.get_ring_for_spin_vec(vec, npts=nphi)
    tht   = np.pi*0.5 - np.arcsin(zs) #latitude
    phi   = np.arctan2(ys, xs)        #longitude
    tpix=hp.ang2pix(nside,tht,phi)
    psmap2.append(l27m24[tpix])


thtmap=[]
phimap=[]
for ring in range(startring,stopring,stepring):
    vec = vecs[ring] 
    xs, ys, zs = qr.maths.get_ring_for_spin_vec(vec, npts=nphi)
    tht   = np.pi*0.5 - np.arcsin(zs) #latitude
    phi   = np.arctan2(ys, xs)        #longitude
    zmin=x[zs==np.min(zs)]
    thtmap.append(roll(tht,zmin[0]))
    phimap.append(roll(phi,zmin[0]))


thtmap2=[]
phimap2=[]
for ring in range(startring,stopring,stepring):
    vec = vecs[ring] 
    xs, ys, zs = qr.maths.get_ring_for_spin_vec(vec, npts=nphi)
    tht   = np.pi*0.5 - np.arcsin(zs) #latitude
    phi   = np.arctan2(ys, xs)        #longitude
    thtpneg=tht[phi<-np.pi/4]
    thtmin=x[np.abs(tht-np.pi/2)==np.min(np.abs(thtpneg-np.pi/2))]
    thtmap2.append(roll(tht,-thtmin[0]))
    phimap2.append(roll(phi,-thtmin[0]))
    
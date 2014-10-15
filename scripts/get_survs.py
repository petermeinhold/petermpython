
from planck.Planck import Planck
import h5py
import healpy as hp
import os
import numpy as np
from glob import glob
"""
script to load surveys and diff 

"""

def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()
    
base='/global/scratch/sd/planck/user/zonca/resdst/'    
freq=44
chan=0    #remember numberin in resdst goes from 0 to n where n is 2*number of radiomters in freq
pol='I'   #choose which maps to get

dirname='dx8s_30'
dirname='dx8s_44'
dirname='dx8s_70'
dirname='dvv_fix_30'
dirname='dvv_fix_44'
dirname='dvv_fix_70'
nside = 1024
#psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
#wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)

#get survey 1
m1= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv1/map.h5',mode='r')['data'][pol]),nside,order_in='NESTED',order_out='RING'))
#get survey 2
m2= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv2/map.h5',mode='r')['data'][pol]),nside,order_in='NESTED',order_out='RING'))
#get survey 3
m3= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv3/map.h5',mode='r')['data'][pol]),nside,order_in='NESTED',order_out='RING'))
#get survey 4
m4= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv4/map.h5',mode='r')['data'][pol]),nside,order_in='NESTED',order_out='RING'))


totalmask=m1.mask|m2.mask|m3.mask|m4.mask#|psmask|wmask

m1.mask=totalmask
m2.mask=totalmask
m3.mask=totalmask
m4.mask=totalmask




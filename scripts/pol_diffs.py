from planck.Planck import Planck
import h5py
import healpy as hp
import os
import numpy as np
from glob import glob

def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()
    
base='/global/scratch/sd/planck/user/zonca/resdst'    
freq=70
dirname='/dx8s_70'
#os.chdir(base+dirname)
    
nside = 1024
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)

#get survey 1-2
q12= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv1-2/map.h5',mode='r')['data']['Q']),nside,order_in='NESTED',order_out='RING'))
u12= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv1-2/map.h5',mode='r')['data']['U']),nside,order_in='NESTED',order_out='RING'))
i12= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv1-2/map.h5',mode='r')['data']['I']),nside,order_in='NESTED',order_out='RING'))

#survyeyt 3-4
q34= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv3-4/map.h5',mode='r')['data']['Q']),nside,order_in='NESTED',order_out='RING'))
u34= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv3-4/map.h5',mode='r')['data']['U']),nside,order_in='NESTED',order_out='RING'))
i34= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+'/surv3-4/map.h5',mode='r')['data']['I']),nside,order_in='NESTED',order_out='RING'))

totalmask=q12.mask|u12.mask|i12.mask|q34.mask|u34.mask|i34.mask|psmask|wmask

q12.mask=totalmask
u12.mask=totalmask
i12.mask=totalmask
q34.mask=totalmask
u34.mask=totalmask
i34.mask=totalmask

#y1=np.reshape(np.concatenate([i12,q12,u12],axis=1),3,len(i12))
#y2=np.reshape(np.concatenate([i34,q34,u34],axis=1),3,len(i12))
#ysum=np.reshape(np.concatenate([(i12+i34)/2.,(q12+q34)/2.,(u12+u34)/2.],axis=1),3,len(i12))
#ydif=np.reshape(np.concatenate([(i12-i34)/2.,(q12-q34)/2.,(u12-u34)/2.],axis=1),3,len(i12))

cl_y1=hp.anafast([i12.filled(),q12.filled(),u12.filled()],nspec=2,lmax=1024)
cl_y2=hp.anafast([i34.filled(),q34.filled(),u34.filled()],nspec=2,lmax=1024)
cl_ysum=hp.anafast([(i12.filled()+i34.filled())/2.,(q12.filled()+q34.filled())/2.,(u12.filled()+u34.filled())/2.],nspec=2,lmax=1024)
cl_ydif=hp.anafast([(i12.filled()-i34.filled())/2.,(q12.filled()-q34.filled())/2.,(u12.filled()-u34.filled())/2.],nspec=2,lmax=1024)



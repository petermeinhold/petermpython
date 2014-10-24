from planck.Planck import Planck
import sys
import os
import numpy as np
import healpy as hp
import planck_util_prm as pu
from glob import glob

if __name__ == "__main__":
    nside =256
    topdir='/global/project/projectdirs/planck/data/ffp8/mc_noise/'
    pkldir='/global/homes/p/peterm/ffp8_noise_cls_1000/'
    tmask=pu.get_lfi_dx11_mask(256)
    mcdirlist=np.sort(glob(topdir+freq+'/'+freq+'_0*'))
    homedir='/global/homes/p/peterm/'
    #arguments: freq s1,q1,hr1,s2, q2, hr2
    freq=sys.argv[1]
    q1=sys.argv[2]
    s1=sys.argv[3]
    hr1=sys.argv[4]
    q2=sys.argv[5]
    s2=sys.argv[6]
    hr2=sys.argv[7]
    
    cls=[]
    for mcdir in mcdirlist:
        #extract the 'hundreds' column to build list of mc numbers in this directory
        h=mcdir[-1]
        nstart=100*int(h)
        mcnumlist=[str(i).zfill(5) for i in range(nstart,nstart+100)]
        for mcnum in mcnumlist:
            f1=mcdir+'/ffp8_noise_%s_%s_map_mc_%s%s.fits' %(freq,s1,mcnum,hr1)
            f2=mcdir+'/ffp8_noise_%s_%s_map_mc_%s%s.fits' %(freq,s2,mcnum,hr2)
            cls.append(pu.read_and_diff_files_fast(f1,f2,nside=256,tmask=tmask))           
    pklfilename=pkldir+'/ffp8_noise_null_cls_'+freq+s1+hr1+s2+hr2+'.pkl'
    pklfile=open(pklfilename,'wb')
    cPickle.dump(cls,pklfile)
    pklfile.close()

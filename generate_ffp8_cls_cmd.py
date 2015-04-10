import sys
import os
import numpy as np
import healpy as hp
from glob import glob
import pickle

def degrade_mask(inmask,nside_out=256):
    #stupid function to degrade a mask by making a map and degrading that
    m=hp.ma(np.ones(len(inmask)))
    m.mask=inmask
    mdg=hp.ud_grade(m,nside_out=nside_out)
    return mdg.mask
    
def get_lfi_dx11_mask(nside,masktype='pol',ps=True):
    """
    now using masks suggested by AZa on 1/27/2015, common mask, should already have PS
    apo=true is apodized, masktype='pol' is polarized mask, masktype='int' intensity mask
    """
    maskdir='/global/homes/p/peterm/masks/'
    f=maskdir+'dx11_v2_common_%s_mask_010a_1024.fits' %masktype
    tmask=hp.ma(hp.read_map(f)) 
    tmask=degrade_mask(tmask,nside_out=nside)
    tmask=logical_not(tmask)
    if ps:
        fpsmask30='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/MASKs/mask_ps_30GHz_beam33amin_nside2048.00_DX9_nonblind_holesize3.fits'
        psmask30 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(fpsmask30), nside)))
        tmask=psmask30.astype(np.bool)|tmask.astype(np.bool)
    return tmask
    
def read_and_diff_files_fast(f1,f2,nside=256,tmask=None,return_map=False):
    #assume tmask input is already degraded
    
    mm1=hp.read_map(f1,[0,1,2],verbose=False)
    mm2=hp.read_map(f2,[0,1,2],verbose=False)

    mmm1=[]
    mmm2=[]
    for m1,m2 in zip(mm1,mm2):
        m1=hp.ud_grade(hp.ma(m1),nside_out=nside)
        m2=hp.ud_grade(hp.ma(m2),nside_out=nside)
        tmask=m1.mask | m2.mask | tmask
        mmm1.append(m1)
        mmm2.append(m2)
    
    diff=[]
    for m1,m2 in zip(mmm1,mmm2):
        d=m1-m2
        d.mask=tmask
        diff.append(d)
    
    skyfrac=1-float(tmask.sum())/len(tmask)
        
    cldata=hp.anafast(diff)
    cldata_out=[]
    for cl in cldata:
        cldata_out.append(cl/skyfrac)
        
    if return_map is False:
        return cldata_out
    if return_map is True:
        return cldata_out,diff

if __name__ == "__main__":
    #arguments: freq s1,q1,hr1,s2, q2, hr2
    freq=sys.argv[1]
    qq1=sys.argv[2]
    s1=sys.argv[3]
    hhr1=sys.argv[4]
    qq2=sys.argv[5]
    s2=sys.argv[6]
    hhr2=sys.argv[7]
    q1=qq1.replace('null','')
    q2=qq2.replace('null','')
    hr1=hhr1.replace('null','')
    hr2=hhr2.replace('null','')
    nside =256
    topdir='/global/project/projectdirs/planck/data/ffp8/mc_noise/'
    pkldir='/global/homes/p/peterm/ffp8_noise_cls_1000_dx11maskps/'
    tmask=get_lfi_dx11_mask(256,masktype='int',ps=True)

    mcdirlist=np.sort(glob(topdir+freq+'/'+freq+'_0*'))
    homedir='/global/homes/p/peterm/'
    
    cls=[]
    for mcdir in mcdirlist:
        #extract the 'hundreds' column to build list of mc numbers in this directory
        h=mcdir[-1]
        nstart=100*int(h)
        mcnumlist=[str(i).zfill(5) for i in range(nstart,nstart+100)]
        for mcnum in mcnumlist:
            f1=mcdir+'/ffp8_noise_%s%s_%s_map_mc_%s%s.fits' %(freq,q1,s1,mcnum,hr1)
            f2=mcdir+'/ffp8_noise_%s%s_%s_map_mc_%s%s.fits' %(freq,q2,s2,mcnum,hr2)
            print(f1)
            cls.append(read_and_diff_files_fast(f1,f2,nside=256,tmask=tmask))           
    pklfilename=pkldir+'/ffp8_noise_null_cls_'+freq+q1+s1+hr1+q2+s2+hr2+'.pkl'
    pklfile=open(pklfilename,'wb')
    pickle.dump(cls,pklfile)
    pklfile.close()

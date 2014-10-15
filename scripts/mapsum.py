import numpy as np
import pycfitsio
import healpy as hp

def mapsum inputfilelist,freq=70,wmapmask=True,psmaskmask=True:
"""function to read in a list of maps and sum them using
a union mask. Default options to include wmap mask and
pt source mask. Assume IQU inputs, read separately, mask
sum and into list default freq(for psmask) =70
"""
    #read one map, just to find out size
    m=hp.ma(hp.read_map(inputfilelist[0]),0)
    nside=hp.npix2nside(len(m))
    tmask=m.mask
    if wmapmask==True:
        wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
        tmask=tmask|wmask
    if psmaskmask==True:
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        tmask=tmask|psmask
    nmaps=len(inputfilelist)
    for f in inputfilelist:
        t=hp.ma(hp.read_map(f,0))
        q=hp.ma(hp.read_map(f,1))
        u=hp.ma(hp.read_map(f,2))
        tmask=tmask|t.mask|q.mask|u.mask
        t.mask=tmask
        q.mask=tmask
        u.mask=tmask
        tsum=tsum+t/nmaps
        qsum=qsum+q/nmaps
        usum=usum+u/nmaps
    outmap=[tsum,qsum,usum]

    return outmap

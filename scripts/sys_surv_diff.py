from glob import glob
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import os
import healpy
from planck import ps
from collections import OrderedDict
import pycfitsio
from planck.LFI import LFI
def degrademask(m, nside):
    return np.ceil(hp.ud_grade(m.astype(np.float), nside)).astype(np.bool)

def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()

fwhm_rad=10.*np.pi/180.0

fl=glob('*.fits')
fbases=[]
for f in fl:
    fbases.append(f[:-13])
fbases=np.unique(fbases)

for f in fbases:
    fss1=f+'survey_1.fits'
    fss2=f+'survey_2.fits'
    ss1=hp.ma(hp.read_map(fss1))
    ss2=hp.ma(hp.read_map(fss2))
    tmask=ss1.mask | ss2.mask
    dif12=(ss1-ss2)/2.0
    dif12.mask=tmask
    dif12sm=hp.smoothing(hp.remove_monopole(dif12),fwhm_rad)
    hp.mollview(dif12sm,title= f+' SS1-SS2 ')
    savefig(f+'SS1-SS2.png')

    
    
NSIDE = 1024
psmask = np.floor(hp.ud_grade(
hp.read_map(glob.glob('/global/homes/z/zonca/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), NSIDE))

gal_cut_deg = 2
gal_cut = np.zeros(hp.nside2npix(NSIDE), dtype=np.bool)
#gal_cut[hp.query_strip(NSIDE, np.radians(90-gal_cut_deg), np.radians(90+gal_cut_deg))] = True

ch = str(f.freq)
#for chnum in range(4): 
chnum = 0
chobj = f.ch[chnum]
#ch = chobj.tag
for chnum,chobj in enumerate(f.ch):  
    print(chobj)
    ch = chobj.tag
    if ch == "LFI23M":
        #for calibration in ["dvv_0_dm", "estsl0_dm"]:
        print "Read maps"
        if calibration.find("DPC")>=0:
            maps = map(pycfitsio.read_map, [glob.glob(os.path.join(dpcbasefolder, 'LFI_%s_1024_*_%s_survey_%d.fits' % (chobj.f.freq, ch.replace('LFI',''), surv)))[0] for surv in survey])
        else:
            if ch.startswith('LFI'):
                maps = []
                for surv in survey:
                    filename = basefolder + '%s/surv%d/map_ch%d.h5' % (calibration,surv, chnum)
                    print(filename)
                    maps.append(hp.ma(hp.ud_grade(maskmap(h5py.File(filename,mode='r')['data'][:]), 1024, order_in='NESTED', order_out='RING')))
            else:
                maps = [hp.ma(hp.ud_grade(h5py.File(basefolder + '%s/surv%d/map.h5' % (calibration,surv ),mode='r')['data']['I'], 1024, order_in='NESTED', order_out='RING')) for surv in survey]

        for m in maps:
            m.mask = m==0
            m -= m.mean()

        print "Survey diff"
        for short, combname, comb in combs:
            print(combname)
            s = comb
            s-=s.mean()
            mask = s.mask | np.logical_not(psmask) | gal_cut
            s.mask = mask
            s = hp.ma(hp.smoothing(s.filled(), fwhm=np.radians(3)))*1e3
            #s = hp.ma(hp.ud_grade(hp.smoothing(s.filled(), fwhm=np.radians(3.)), 1024))*1e6
            s.mask = mask   
            hp.mollview(s.filled(),min=-.1,max=.1,unit='mK', title='%s: %s %s' % (ch, combname, calibration))
            plt.savefig('pngsurvdiff/' + '_'.join([calibration, short, ch]) + '.png')

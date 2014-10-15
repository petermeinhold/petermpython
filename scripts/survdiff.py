import glob
import h5py
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import itertools
import os
import healpy
from planck import ps
from collections import OrderedDict
import pycfitsio
from planck.LFI import LFI

lfi = LFI()

freq = 70
f = lfi.f[freq]

dpcbasefolder = "/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi"
dpcbasefolder = "/global/scratch/sd/planck/user/zonca/DX8_dvv/"
basefolder = "/global/scratch/sd/planck/user/zonca/resdst/"

survey = [1,2,3]
survey = [1]#,2]
survey = [1,2]
survey = [1,2,3,4]

calibration = 'dvv_tom_30'
calibration = 'hk_30'
calibration = 'dvv_suf_30_td'
calibration = 'predx9s_30_td'
calibration = 'dvv_tom_30_td'
calibration = 'DPC_dvv'
calibration = 'dvv_tfix_%d_td' % freq
calibration = 'dx8s_%d_td' % freq
calibration = 'dx8s_%d_nodm_td' % freq
calibration = 'dvv_fix_%d_nor_td' % freq
#for calibration in ["dstcalsmt_dm_ip", "DX8S_dm_ip"]:

def degrademask(m, nside):
    return np.ceil(hp.ud_grade(m.astype(np.float), nside)).astype(np.bool)

def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()

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
        #diffmaps = {}
        #for i,j in itertools.combinations(survey,2):
        #    print(i,j)
        #    if (i,j) == (2, 3):
        #        i, j = 3, 2
        #    print(i,j)
        #    diff = maps[survey.index(j)]- maps[survey.index(i)] 
        #    diff -= diff.mean()
        #    diffm = ps.smooth(diff.filled(), 60.)
        #    diffmaps[(i,j)] = hp.ma(diffm)*1e6
        #    hp.mollview((diffmaps[(i,j)]).filled(), min=-50,max=50,unit='muK', title='%s: SS%d - SS%d %s' % (ch, j, i, calibration))
        #    plt.savefig('png/%s_%s_SS%d_SS%s.png' % (calibration, ch, j, i))
        combs =( #('SS24-SS13', '(SS4+SS2)/2 - (SS1+SS3)/2',((maps[1]+maps[3])-(maps[2]+maps[0]))/2.),
                 ('SS2-SS1', '(SS2-SS1)',maps[1]-maps[0]),
                 ('SS4-SS3', '(SS4-SS3)',maps[3]-maps[2]),
                 ('SS4-SS2', '(SS4-SS2)',maps[3]-maps[1]),
                 ('SS4-SS1', '(SS4-SS1)',maps[3]-maps[0]),
                 #('SS3-SS1', '(SS3-SS1)',maps[2]-maps[0]),
                 #('SS3-SS2', '(SS3-SS2)',maps[2]-maps[1]),
                 )

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
            #np.save('pklsurvdiff/' + '_'.join([calibration, short, ch]) + '.npy', s.filled())
            #dst_map = hp.ma(hp.ud_grade(s.filled(), 1024, order_in='RING', order_out='NESTED'))/1.e6
            #dst_map += np.abs(dst_map.min())
            #dst_map.fill_value = 0 # for file
            #with h5py.File('_'.join([calibration, short, ch]) + '.h5', mode='w') as h5file:
            #    h5file.create_dataset("data", data=dst_map.filled())

        #print("Sidelobes")
        #fslfolder = '/global/homes/z/zonca/gpeterm/sidelobes/perrotta/'
        ##fslmaps = [hp.ma(hp.read_map(fslfolder + 'madam_sidelobes_galaxy_30GHz_ss%d_512binmap.fits' % surv)) for surv in [1,2]]
        #fslmaps = [hp.ma(hp.read_map(fslfolder + 'madam_sidelobes_galaxy_30GHz_ss%d_512outmap.fits' % surv)) for surv in [1,2]]
        #fsldiff = fslmaps[1]-fslmaps[0]
        #fsldiff -= fsldiff.mean()
        #fsldiff.mask = fsldiff.mask | degrademask(s.mask, 512)
        #fsldiffm = fsldiff*1e6
        #fsldiffm = hp.ma(ps.smooth(fsldiff.filled(), 60.))*1e6
        #hp.mollview(fsldiffm.filled(),min=-12,max=12,unit='muK', title='FSL SIM (Perrotta) %s: (SS2-SS1) %s' % (ch, calibration))
        #plt.text(1.25,-.9,'Zonca - Meinhold')
        #plt.savefig('fsl.png')
        #hp.mollview((s-fsldiffm/12.*50).filled(),min=-50,max=50,unit='muK', title='Data - scaledSIM')
        #plt.text(1.25,-.9,'Zonca - Meinhold')
        #plt.savefig('res.png')

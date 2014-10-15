import numpy as np
from planck import ps
import os
import pycfitsio
import glob
import healpy as hp
import matplotlib.pyplot as plt
import h5py

NSIDE = 16
mask = hp.ud_grade(hp.read_map('../singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),NSIDE) < 1
basefolder = "/global/scratch/sd/planck/user/zonca/resdst/"

calibration = "DX8bR_smt"
calibration = "dstcalsmt_dstmask"
calibration = "DX8S"
calibration = "dstcalsmt_dm_ip"
calibration = 'cal_dm'
calibration = 'DX8S_dm'
calibration = "dvv_mf_dm_fsl"
calibration = 'dvv_tom_dm'
calibration = "flat_dm"
calibration = "DX8S_30_fsl"
calibration = 'dvv_tom_fsl'
calibration = "dvv_mf_dm"
calibration = "DX8S_30"
calibration = 'dvv_mf_30'
calibration = 'hk_30'
calibration = 'predx9s_30'
calibration = 'dvv_tom_30'
calibration = 'predx9s_fix_30'
calibration = 'dvv_sky_30'
#for freq in [30,44,70]:#freq = 70
cor_folder = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/qucorrection_dx8/v1/"
cor = False
fit_factor = False

central_freqs = { 30:29.5, 44:44.1, 70:70.3 }
WMAP_freqs = { 'K':22.8, 'Ka':33., 'Q':40.7, 'V':60.8, 'W':93.5 }
sync_spectral_index = -3.1

for freq in [30]:#freq = 70
    calibration = 'dvv_fix_%d' % freq
    calibration = "DPCdvv_%d" % freq
    calibration = "DPCDX8_%d" % freq

    if calibration.find("DX8")>0:
        dpcbasefolder = "/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi"
    else:
        dpcbasefolder = "/project/projectdirs/planck/data/mission/DPC_maps/predx9/lfi"
    if cor:
        calibration = calibration.replace("_%d" % freq, "_cor_%d" % freq)

    def map_fit(m1, m2):
        return np.sum(m1.filled() * np.logical_not(m1.mask|m2.mask) \
        * m2.filled())/np.sum((m2 ** 2).filled() * np.logical_not(m1.mask|m2.mask))
    ch = freq
    it = 39
    it = 0
    #for surv in [0]:#1,2]:#,3,4]:
    #for surv in [0]:
    #for surv in [1,2,3,4]:
    for num,comp in enumerate("Q"):
        k=hp.ma(hp.ud_grade(hp.read_map('../singler/wmap/wmap_band_iqumap_r9_7yr_K_v4.fits',num+1), NSIDE))*1e3
        for surv in [0]:

            if calibration.find("DPC")>=0:
                filename = 'LFI_%s_1024_????????.fits' % ch
                if surv > 0:
                    filename = filename.replace(".fits", "_survey_%d.fits" % surv)
                else:
                    if calibration.find("dvv")>0:
                        filename = filename.replace(".fits", "_dx8.fits")
                filename_found = glob.glob(os.path.join(dpcbasefolder, filename))[0]
                print filename_found
                qmap = hp.read_map(filename_found,num+1) 
                if cor:
                    qmap += hp.read_map(cor_folder + "iqu_bandpass_correction_%d_fullsurvey.fits" % freq, num + 1)
                m = hp.ma(hp.ud_grade(qmap, NSIDE))*1e6
            else:
                #qmap = hp.ma(h5py.File(basefolder + '%s/surv%d/%02d_map.h5' % (calibration,surv, it))['data']['Q'])
                qmap = hp.ma(h5py.File(basefolder + '%s/surv%d/map.h5' % (calibration,surv))['data'][comp])
                qmap.mask = qmap == 0
                if cor:
                    qmap += hp.ma(hp.read_map(cor_folder + "iqu_bandpass_correction_%d_fullsurvey.fits" % freq, num + 1))
                qmap = qmap.filled()
                m = hp.ma(hp.ud_grade(qmap, NSIDE, order_in='NESTED', order_out='RING',pess=False))*1e6

                #radius = np.radians(chobj.get_instrument_db_field('beamfwhm') * 2)
                #total_mask[
                #healpy.query_disc(
                #   1024, healpy.ang2vec( np.radians(90-62.4), np.radians(-88.2) ), radius 
                #   )] = True

            k.mask = mask | k.mask
            m.mask = mask | m.mask | (m==0)

            fac = (central_freqs[freq]/WMAP_freqs['K'])**sync_spectral_index
            fac = .5
            if fit_factor:
                fac = map_fit(m,k)
            print(fac)
            diff = (m-fac*k)
            #diffm = ps.smooth(diff.filled(), 30.)
            #diffm = hp.ud_grade(diff.filled(), 16)
            diffm = diff.filled()
            hp.mollview(diffm, title=('%s ' + comp + ' - %.4g Kband '+comp+', Survey %d') % (calibration, fac, surv),min=-40,max=40,  unit='%')
            hp.visufunc.projtext(np.pi/2,np.radians(100),'PeaktoPeak %.2f muK, RMS %.2f muK' % (np.ptp(diff), np.std(diff)),size='x-large')
            plt.savefig('polwmap/%s_%s_%d.png' % (calibration, comp,surv))

    #plt.figure()
    #plt.hist(dmr.compressed(), bins=50, range=[-30,30], log=True, label='MultiR')
    #plt.hist(dsr.compressed(), bins=50, range=[-30,30], log=True, label='SingleR')
    #plt.xlabel('muK'); plt.ylabel('Hitcount'); plt.title('Difference map histogram')
    #plt.legend()
    #plt.savefig('histogram.png')

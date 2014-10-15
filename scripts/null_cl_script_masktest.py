import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu
import cPickle
basedelta='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'

#explicit filenames to avoid confusion
f30rh1= '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_1_nominal_1s.fits'
f30rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_2_nominal_1s.fits'
f30ss1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_survey_1_1s.fits'
f30ss2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_survey_2_1s.fits'
f30ss1rh1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_1_survey_1_1s.fits'
f30ss1rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_2_survey_1_1s.fits'
f30ss2rh1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_1_survey_2_1s.fits'
f30ss2rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_2_survey_2_1s.fits'

nside=1024
psmask30 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_30GHz_*.*')[0]), nside)))

mask80 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/galactic_new/*.80_sky*.fits')[0]), nside)))

m30ss1=hp.ma(hp.read_map(f30ss1,[0,1,2])) 
m30ss2=hp.ma(hp.read_map(f30ss2,[0,1,2]))

m30ssdiff=[]
tmask30=psmask30.astype(np.bool)|mask80.astype(np.bool)
m1mask=copy(m30ss1[1].mask)
m2mask=copy(m30ss2[1].mask)
for m1,m2 in zip (m30ss1,m30ss2):
    m1.mask=tmask30|m1mask|m2mask
    m2.mask=tmask30|m1mask|m2mask
    m30ssdiff.append(m1-m2)
    frac30ss=(1-float(m1.mask.sum())/len(m1.mask))     

cl30ssdiff_ptsrc_gal=np.array(hp.anafast(m30ssdiff))*1e12/frac30ss
m30ssdiff=[]
for m1,m2 in zip (m30ss1,m30ss2):
    m1.mask=psmask30|m1a.mask|m2a.mask
    m2.mask=psmask30|m1a.mask|m2a.mask
    m30ssdiff.append(m1-m2)
    frac30ss=(1-float(m1.mask.sum())/len(m1.mask))     

cl30ssdiff_ptsrc=np.array(hp.anafast(m30ssdiff))*1e12/frac30ss
m30ssdiff=[]
for m1,m2 in zip (m30ss1,m30ss2):
    m1.mask=mask80|m1a.mask|m2a.mask
    m2.mask=mask80|m1a.mask|m2a.mask
    m30ssdiff.append(m1-m2)
    frac30ss=(1-float(m1.mask.sum())/len(m1.mask))     

cl30ssdiff_gal=np.array(hp.anafast(m30ssdiff))*1e12/frac30ss

m30ssdiff=[]
for m1,m2 in zip (m30ss1,m30ss2):
    m1.mask=m1a.mask|m2a.mask
    m2.mask=m1a.mask|m2a.mask
    m30ssdiff.append(m1-m2)
    frac30ss=(1-float(m1.mask.sum())/len(m1.mask))     

cl30ssdiff_none=np.array(hp.anafast(m30ssdiff))*1e12/frac30ss

l=np.arange(len(cl30ssdiff_ptsrc_gal[0]))
figure()
plot(cl30ssdiff_ptsrc_gal[0],label='ssdiff, ptsrc and gal masks' )
plot(cl30ssdiff_ptsrc[0],label='ssdiff, ptsrc mask' )
plot(cl30ssdiff_gal[0],label='ssdiff, gal mask' )
plot(cl30ssdiff_none[0],label='ssdiff, no extra mask' )
xlabel('Multipole'),ylabel('Cl, corrected for skyfraction')
xscale('log'),yscale('log')
title='Masking of null test survey differences SS1-SS2 30 GHz'
leg=plt.legend()
pu.thicklegendlines(leg)
xlim([1,3000])
ylim([.01,100])
grid()
plt.savefig('/global/homes/p/peterm/cl_null_30_masks.png')


import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu
import cPickle
basedelta='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
#UPDATED MARCH 18, 2013 to include afactor of 1/2 for all diff maps, to match mean maps (as for null tests)

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
psmask44 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_44GHz_*.*')[0]), nside)))
psmask70 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_70GHz_*.*')[0]), nside)))

mask80 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/galactic_new/*.80_sky*.fits')[0]), nside)))
tmask30=psmask30.astype(np.bool)|mask80.astype(np.bool)

m30ss1=hp.ma(hp.read_map(f30ss1,[0,1,2])) 
m30ss2=hp.ma(hp.read_map(f30ss2,[0,1,2]))
#m30ss1_512=hp.ud_grade(m30ss1,nside_out=nside)
#m30ss2_512=hp.ud_grade(m30ss2,nside_out=nside)
#m30ssdiff_512=[]
m30ssdiff=[]
for m1,m2 in zip (m30ss1,m30ss2):
    m1.mask=tmask30|m1.mask|m2.mask
    m2.mask=tmask30|m1.mask|m2.mask
    m30ssdiff.append((m1-m2)/2.)
    frac30ss=(1-float(m1.mask.sum())/len(m1.mask))     
cl30ssdiff=np.array(hp.anafast(m30ssdiff))*1e12/frac30ss
l=np.arange(len(cl30ssdiff[0]))


figure()
plot l,cl30ssdiff[0]*l*(l+1)/(2*np.pi),label='PM ssdiff'
plot aml,1e12*amcl*aml*(aml+1)/(2*np.pi),label='AM ssdiff'
plot aml,1e12*amcl*aml*(aml+1)/(2*np.pi)-cl30ssdiff[0][:1536]*aml*(aml+1)/(2*np.pi),label='AM-PM ssdiff'
leg=legend()
grid()
xlim([0,800])
ylim([-100,1000])
xlabel('Multipole',fontsize=18)
ylabel('$l(l+1)C_l/2 \pi , \mu K^2$',fontsize=18)
title('30 GHz (SS1-SS2)/2 comparisons',fontsize=20)
show()

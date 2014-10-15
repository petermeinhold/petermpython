import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util as pu
import cPickle

#explicit filenames to avoid confusion

f30ss1='/global/project/projectdirs/planck/data/releases/LFI_PR1/LFI_SkyMap_030_1024_R1.10_survey_1.fits'
f30ss2='/global/project/projectdirs/planck/data/releases/LFI_PR1/LFI_SkyMap_030_1024_R1.10_survey_2.fits'
fpsmask30='/global/project/projectdirs/planck/data/releases/LFI_PR1/LFI_MASK_030-ps_2048_R1.00.fits'
fmask70='/global/project/projectdirs/planck/data/releases/LFI_PR1/COM_MASK_gal-07_2048_R1.00.fits'

nside=1024
#read and degrade the gal an pt src masks

psmask30 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(fpsmask30), nside)))
mask70 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(fmask70), nside)))

#combine the galaxy and pt src masks
tmask30=psmask30.astype(np.bool)|mask70.astype(np.bool)

#read the maps, make them masked arrays
m30ss1=hp.ma(hp.read_map(f30ss1)) 
m30ss2=hp.ma(hp.read_map(f30ss2))

#union of al the masks
m30ss1.mask=tmask30|m30ss1.mask|m30ss2.mask
m30ss2.mask=tmask30|m30ss1.mask|m30ss2.mask
m30ssdiff=(m30ss1-m30ss2)/2.

#output some maps, just to verify the masks and overall agreement
hp.mollview m30ss1,min=-1e-3,max=1e-3,title='PR1 30 GHz survey 1 70% gal mask'
savefig('plots/clplots/pr_map_30_ss1.png')
hp.mollview m30ss1,min=-1e-3,max=1e-3,title='PR1 30 GHz survey 2 70% gal mask'
savefig('plots/clplots/pr_map_30_ss2.png')
hp.mollview m30ss1,min=-1e-4,max=1e-4,title='PR1 30 GHz (SS1 - SS2)/2 70% gal mask'
savefig('plots/clplots/pr_map_30_ss12diff.png')

#calculate and print sky fraction
frac30ss=(1-float(m30ss1.mask.sum())/len(m30ss1.mask))
print frac30ss

#calculate spectrum, adjust for sky fraction and to microK^2
cl30ssdiff=hp.anafast(m30ssdiff)*1.e12/frac30ss

#dump the spectrum to pickle file
fcl=open('plots/clplots/pr1_cl_30ghzssdiff.pkl','wb')
cPickle.dump(cl30ssdiff,fcl)

#plot the spectrum
l=np.arange(len(cl30ssdiff)) + 1
figure()
plot l,cl30ssdiff*l*(l+1)/(2*np.pi),label='PM ssdiff'
leg=legend()
grid()
xlim([0,800])
ylim([-100,1000])
yscale('log')
xlabel('Multipole',fontsize=18)
ylabel('$l(l+1)C_l/2 \pi , \mu K^2$',fontsize=18)
title('Public Release 1 30 GHz (SS1-SS2)/2 comparisons',fontsize=20)
show()
savefig('plots/clplots/pr1_cl_30ghzssdiff.png')


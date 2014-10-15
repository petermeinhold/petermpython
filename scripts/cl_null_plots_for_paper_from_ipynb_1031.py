# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
from petermpython import planck_util_prm as pu
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

f44rh1= '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_1_nominal_1s.fits'
f44rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_2_nominal_1s.fits'
f44ss1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_survey_1_1s.fits'
f44ss2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_survey_2_1s.fits'
f44ss1rh1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_1_survey_1_1s.fits'
f44ss1rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_2_survey_1_1s.fits'
f44ss2rh1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_1_survey_2_1s.fits'
f44ss2rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_2_survey_2_1s.fits'

f70rh1 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_1_nominal_1s.fits'
f70rh2 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_2_nominal_1s.fits'
f70ss1 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_survey_1_1s.fits'
f70ss2 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_survey_2_1s.fits'
f70ss1rh1 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_ringhalf_1_survey_1_1s.fits'
f70ss1rh2 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_ringhalf_2_survey_1_1s.fits'
f70ss2rh1 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_ringhalf_1_survey_2_1s.fits'
f70ss2rh2 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_ringhalf_2_survey_2_1s.fits'

nside=1024

# <codecell>

psmask30 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_30GHz_*.*')[0]), nside)))
psmask44 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_44GHz_*.*')[0]), nside)))
psmask70 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_70GHz_*.*')[0]), nside)))

mask80 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/galactic_new/*.80_sky*.fits')[0]), nside)))
tmask30=psmask30.astype(np.bool)|mask80.astype(np.bool)
tmask44=psmask44.astype(np.bool)|mask80.astype(np.bool)
tmask70=psmask70.astype(np.bool)|mask80.astype(np.bool)

m30rh1=hp.ma(hp.read_map(f30rh1,[0,1,2]))
m30rh2=hp.ma(hp.read_map(f30rh2,[0,1,2]))
m30ss1=hp.ma(hp.read_map(f30ss1,[0,1,2])) 
m30ss2=hp.ma(hp.read_map(f30ss2,[0,1,2]))

m30ss1rh1=hp.ma(hp.read_map(f30ss1rh1,[0,1,2])) 
m30ss1rh2=hp.ma(hp.read_map(f30ss1rh2,[0,1,2]))
m30ss2rh1=hp.ma(hp.read_map(f30ss2rh1,[0,1,2])) 
m30ss2rh2=hp.ma(hp.read_map(f30ss2rh2,[0,1,2]))

m30rhdiff=[]
m30ssdiff=[]
m30ssrhdiff=[]
for m1,m2 in zip(m30rh1,m30rh2):
    m1.mask=tmask30|m1.mask|m2.mask
    m2.mask=tmask30|m1.mask|m2.mask
    m30rhdiff.append(hp.remove_monopole(m1-m2)/2.)
    frac30rh =(1-float(m1.mask.sum())/len(m1.mask))     
for m1,m2 in zip (m30ss1,m30ss2):
    m1.mask=tmask30|m1.mask|m2.mask
    m2.mask=tmask30|m1.mask|m2.mask
    m30ssdiff.append(hp.remove_monopole(m1-m2)/2.)
    frac30ss=(1-float(m1.mask.sum())/len(m1.mask))     
for m1a,m1b,m2a,m2b in zip (m30ss1rh1,m30ss1rh2,m30ss2rh1,m30ss2rh2):
    m1a.mask=tmask30|m1a.mask|m2a.mask|m1b.mask|m2b.mask
    m1b.mask=m1a.mask
    m2a.mask=m1a.mask
    m2b.mask=m1a.mask
    m30ssrhdiff.append(hp.remove_monopole(m1a+m2a-m1b-m2b)/4.)
    frac30ssrh=(1-float(m1a.mask.sum())/len(m1a.mask)) 

# <codecell>

cl30ssrhdiff=np.array(hp.anafast(m30ssrhdiff,regression=False))*1e12/frac30ss
cl30rhdiff=np.array(hp.anafast(m30rhdiff,regression=False))*1e12/frac30rh
cl30ssdiff=np.array(hp.anafast(m30ssdiff,regression=False))*1e12/frac30ss
cl30rhx=np.array(hp.anafast(m30rh1,m30rh2,regression=False))*1e12/frac30rh  
l=np.arange(len(cl30ssdiff[0]))


m44rh1=hp.ma(hp.read_map(f44rh1,[0,1,2]))
m44rh2=hp.ma(hp.read_map(f44rh2,[0,1,2]))
m44ss1=hp.ma(hp.read_map(f44ss1,[0,1,2])) 
m44ss2=hp.ma(hp.read_map(f44ss2,[0,1,2]))

m44ss1rh1=hp.ma(hp.read_map(f44ss1rh1,[0,1,2])) 
m44ss1rh2=hp.ma(hp.read_map(f44ss1rh2,[0,1,2]))
m44ss2rh1=hp.ma(hp.read_map(f44ss2rh1,[0,1,2])) 
m44ss2rh2=hp.ma(hp.read_map(f44ss2rh2,[0,1,2]))

m44rhdiff=[]
m44ssdiff=[]
m44ssrhdiff=[]
for m1,m2 in zip(m44rh1,m44rh2):
    m1.mask=tmask44|m1.mask|m2.mask
    m2.mask=tmask44|m1.mask|m2.mask
    m44rhdiff.append(hp.remove_monopole(m1-m2)/2.)
    frac44rh =(1-float(m1.mask.sum())/len(m1.mask))     
for m1,m2 in zip (m44ss1,m44ss2):
    m1.mask=tmask44|m1.mask|m2.mask
    m2.mask=tmask44|m1.mask|m2.mask
    m44ssdiff.append(hp.remove_monopole(m1-m2)/2.)
    frac44ss=(1-float(m1.mask.sum())/len(m1.mask))     
for m1a,m1b,m2a,m2b in zip (m44ss1rh1,m44ss1rh2,m44ss2rh1,m44ss2rh2):
    m1a.mask=tmask44|m1a.mask|m2a.mask|m1b.mask|m2b.mask
    m1b.mask=m1a.mask
    m2a.mask=m1a.mask
    m2b.mask=m1a.mask
    m44ssrhdiff.append(hp.remove_monopole(m1a+m2a-m1b-m2b)/4.)
    frac44ssrh=(1-float(m1a.mask.sum())/len(m1a.mask))     

cl44ssrhdiff=np.array(hp.anafast(m44ssrhdiff,regression=False))*1e12/frac44ss
cl44rhdiff=np.array(hp.anafast(m44rhdiff,regression=False))*1e12/frac44rh
cl44ssdiff=np.array(hp.anafast(m44ssdiff,regression=False))*1e12/frac44ss
cl44rhx=np.array(hp.anafast(m44rh1,m44rh2,regression=False))*1e12/frac44rh

m70rh1=hp.ma(hp.read_map(f70rh1,[0,1,2]))
m70rh2=hp.ma(hp.read_map(f70rh2,[0,1,2]))
m70ss1=hp.ma(hp.read_map(f70ss1,[0,1,2])) 
m70ss2=hp.ma(hp.read_map(f70ss2,[0,1,2]))

m70ss1rh1=hp.ma(hp.read_map(f70ss1rh1,[0,1,2])) 
m70ss1rh2=hp.ma(hp.read_map(f70ss1rh2,[0,1,2]))
m70ss2rh1=hp.ma(hp.read_map(f70ss2rh1,[0,1,2])) 
m70ss2rh2=hp.ma(hp.read_map(f70ss2rh2,[0,1,2]))

m70rhdiff=[]
m70ssdiff=[]
m70ssrhdiff=[]
for m1,m2 in zip(m70rh1,m70rh2):
    m1.mask=tmask70|m1.mask|m2.mask
    m2.mask=tmask70|m1.mask|m2.mask
    m70rhdiff.append(hp.remove_monopole(m1-m2)/2.)
    frac70rh =(1-float(m1.mask.sum())/len(m1.mask))     
for m1,m2 in zip (m70ss1,m70ss2):
    m1.mask=tmask70|m1.mask|m2.mask
    m2.mask=tmask70|m1.mask|m2.mask
    m70ssdiff.append(hp.remove_monopole(m1-m2)/2.)
    frac70ss=(1-float(m1.mask.sum())/len(m1.mask))     
for m1a,m1b,m2a,m2b in zip (m70ss1rh1,m70ss1rh2,m70ss2rh1,m70ss2rh2):
    m1a.mask=tmask70|m1a.mask|m2a.mask|m1b.mask|m2b.mask
    m1b.mask=m1a.mask
    m2a.mask=m1a.mask
    m2b.mask=m1a.mask
    m70ssrhdiff.append(hp.remove_monopole(m1a+m2a-m1b-m2b)/4.)
    frac70ssrh=(1-float(m1a.mask.sum())/len(m1a.mask))     

cl70ssrhdiff=np.array(hp.anafast(m70ssrhdiff,regression=False))*1e12/frac44ss
cl70rhdiff=np.array(hp.anafast(m70rhdiff,regression=False))*1e12/frac44rh
cl70ssdiff=np.array(hp.anafast(m70ssdiff,regression=False))*1e12/frac44ss
cl70rhx=np.array(hp.anafast(m70rh1,m70rh2,regression=False))*1e12/frac44rh


fcl=open('/global/homes/p/peterm/plots/clplots/cl_nulls_nomono.pkl','wb')
cl_nulls={'cl30ssrhdiff':cl30ssrhdiff,'cl30rhdiff':cl30rhdiff,'cl30rhx':cl30rhx,
'cl30ssdiff':cl30ssdiff,'cl44ssrhdiff':cl44ssrhdiff,'cl44rhdiff':cl44rhdiff,'cl44rhx':cl44rhx,
'cl44ssdiff':cl44ssdiff,'cl70ssrhdiff':cl70ssrhdiff,'cl70rhdiff':cl70rhdiff,'cl70rhx':cl70rhx,
'cl70ssdiff':cl70ssdiff}
cPickle.dump(cl_nulls,fcl)
fcl.close()

# <codecell>

import sys
sys.path.append('/global/homes/p/peterm/petermpython/paperplots/python/scripts')
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import cPickle
import numpy as np
import pyfits

# <codecell>

cd ~/plots/clplots

# <codecell>

# Load data
power_spectrum= np.loadtxt("bf_cbipap5_all.dat") # l, C_l
# load WF filtered data (from Andrea Zonca)

f30=pyfits.open('planck_tauprior_highL_post_lensing_beam_30.fits')
f44=pyfits.open('planck_tauprior_highL_post_lensing_beam_44.fits')
f70=pyfits.open('planck_tauprior_highL_post_lensing_beam_70.fits')
ps30=f30[1].data['TEMPERATURE']
ps44=f44[1].data['TEMPERATURE']
ps70=f70[1].data['TEMPERATURE']


fcl=open('cl_nulls_nomono.pkl','rb')
cl_nulls=cPickle.load(fcl)
fcl.close()
l=np.arange(len(cl_nulls['cl30ssdiff'][0]))
# Configure Matplotlib options
# Create the plot
fnames={}
fnames['30']=['cl30rhx','cl30ssdiff','cl30ssrhdiff','cl30rhdiff']
fnames['44']=['cl44rhx','cl44ssdiff','cl44ssrhdiff','cl44rhdiff']
fnames['70']=['cl70rhx','cl70ssdiff','cl70ssrhdiff','cl70rhdiff']

lab=['Ring Half cross spectrum','Survey 1-2 Difference','Ring Half difference, survey 1+2','Ring Half difference, nominal']
logplot=True

# <codecell>

type(ps30)

# <codecell>

for freq,psfreq in zip(['30','44','70'],[ps30,ps44,ps70]):
    names=fnames[freq]    
    for width in [8.8]:#[18., 12., 8.8]:
        fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
        # this should be changed for making a panel of multiple figures
        ax = fig.add_subplot(111)
           
        # power spectrum
        #plt.plot(power_spectrum[:,0], power_spectrum[:,1], "k", label="Best fit")
        plt.plot(psfreq, "k", label="Best fit x window function",alpha=.5)

        for i in [1,2,3]:
            plt.plot(l,l*(l+1)*cl_nulls[names[i]][0]/(2.*np.pi),label=lab[i])
            
        plt.plot(l,l*(l+1)*np.abs(cl_nulls[names[1]][0]-cl_nulls[names[2]][0])/(2.*np.pi),label='Survey 1-2 signal residual')
        
        # legend
        plt.legend(frameon=False,loc=0)

        # labels
        plt.xlabel(r"$\ell$"); plt.ylabel(r"$\ell(\ell+1)C_\ell / 2\pi \:  \left[ \mu K^2 \right]$")
        ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

        # reduce ticks for small figures
        if width < 10:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        
        # axes limits
        if logplot==True:
            plt.ylim([1, 6000]); plt.xlim([10, 800]);
            plt.xscale('log')
            plt.yscale('log')
            # legend
            plt.legend(frameon=False,loc=2)
        if logplot==False:
            plt.ylim([0, 6000]); plt.xlim([0, 800]);
            plt.xscale('linear')
            plt.yscale('linear')
            # legend
            plt.legend(frameon=False,loc=2)
        # reduce white space around figure
        plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

        # set vertical y axis ticklables
        for ticklabel in ax.yaxis.get_ticklabels():
            ticklabel.set_rotation("vertical")
        if logplot==False:
            plt.savefig("null_cls_nomono_"+freq+"ghz_wf%dmm.pdf" % int(width*10), bbox_inches='tight')
        if logplot==True:
            plt.savefig("null_cls_nomono_"+freq+"ghz_log_wf%dmm.pdf" % int(width*10), bbox_inches='tight')

        # <codecell>

# <codecell>



#gain from noise and tys    
# script to retrieve pkl ed tsys and uncalibrated noise estimates (per PID) to generate gain estimates

import planck_util as pu
from glob import glob
import pyfits
import cPickle
from scipy.stats import nanmean,nanstd,nanmedian
from testenv.todtools import write_calfile
import numpy as np
import matplotlib.pyplot as plt

#set up the output moving average:
binsize=30
w = np.empty((binsize,), dtype=np.float_)
w[:] = 1.0/binsize

ucds30=pu.get_ucds(freq=30)
ucds44=pu.get_ucds(freq=44)
ucds70=pu.get_ucds(freq=70)
ucds=dict(ucds70,**ucds44)
ucds.update(ucds30)

ftsys=open('tsys_y_all_fixed.pkl','rb')
tsys_all=cPickle.load(ftsys)
ftsys.close()
def removeoutliers(inarray,stdcut=3.0):
    #bonehead outlier cut, stdcut is how many sigma, replace with nearest neighbor
    #first mark the bad numbers
    inarray[np.logical_not(np.isfinite(inarray))]=0.
    indexarray=np.arange(len(inarray))
    badi=indexarray[np.abs(inarray-nanmedian(inarray)) > stdcut*nanstd(inarray) ]
    goodi=indexarray[np.abs(inarray-nanmedian(inarray)) <= stdcut*nanstd(inarray) ]
    outarray=inarray
    for i in badi:
        outarray[i]=inarray[np.abs(goodi-i).argmin()]
    return outarray

rcalist=[]
for r in range(18,29):
    for d in ['M','S']:
        rcalist.append('LFI'+str(r)+d)
        
        
#70 GHz first
freq=70
tsky=.51*2.735
cols=[]
chcalib=[]
choff=[]
plt.figure(figsize=[12,5])
#need to find universal intersection PIDs among 70 Ghz radiometers:
pid70=ucds['LFI18M']['pID']
for rca in rcalist[:13]:
    pid70=np.intersect1d(pid70,ucds[rca]['pID'])

for rca in rcalist[:12]:
    #get average nominal gain to match
    mediannominalgain=nanmedian(ucds[rca]['nominal_gain'])
    #deparse to the uncal noise file name
    fnoise=open('wn_'+rca+'.pkl','rb')
    wnoise=cPickle.load(fnoise)
    fnoise.close()
    pidout=[]
    noisegainout=[]
    for pid in pid70:
        wn=wnoise['wn7'][np.where(wnoise['pID']==pid)[0]]
        noisegainunscaled=(tsys_all[rca]['tsys_y'][tsys_all[rca]['pID']==pid]+tsky)/wn
        #require calibrated noise to match Tsys*constant (=1/sqrt(bt))
        if len(noisegainunscaled)>0:
            pidout.append(pid)
            noisegainout.append(noisegainunscaled)            
    pidout=np.array(pidout)
    noisegainout=np.array(noisegainout)
    scale=mediannominalgain/nanmedian(noisegainout)
    noisegainout=noisegainout*scale
    noisegainout=removeoutliers(noisegainout.flatten())
    smoothednoisegainout=np.convolve(noisegainout, w, mode='same')
    smoothednoisegainout[:binsize/2]=smoothednoisegainout[binsize/2+1]
    smoothednoisegainout[-binsize/2:]=smoothednoisegainout[-binsize/2-1]
    #make a plot and save it
    plt.hold(False)
    plt.plot(pidout,noisegainout,label='Noise gain')
    plt.hold(True)
    plt.plot(pidout,smoothednoisegainout,label='Binned '+str(binsize),color='black')
    plt.plot(ucds[rca]['nominal_gain'],label='UCDS nominal gain',color='red',lw=2)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.xlabel('Pointing  ID'),plt.ylabel('Gain, K/V'),plt.title(rca + ' Noise estimated gain using Y factor Tsys')
    plt.show()
    plt.savefig('noisegain_'+rca+'.png')
    
    chcalib.append(smoothednoisegainout)
    choff.append(np.zeros(len(smoothednoisegainout)))
write_calfile('tswn',70,pid70,rcalist[:13],chcalib,choff)


#44 GHz next
freq=44
tsky=.66*2.735
cols=[]
chcalib=[]
choff=[]

#need to find universal intersection PIDs among 70 Ghz radiometers:
pid44=ucds['LFI24M']['pID']
for rca in rcalist[13:19]:
    pid44=np.intersect1d(pid44,ucds[rca]['pID'])

for rca in rcalist[12:18]:
    #get average nominal gain to match
    mediannominalgain=nanmedian(ucds[rca]['nominal_gain'])
    #deparse to the uncal noise file name
    fnoise=open('wn_'+rca+'.pkl','rb')
    wnoise=cPickle.load(fnoise)
    fnoise.close()
    pidout=[]
    noisegainout=[]
    for pid in pid44:
        wn=wnoise['wn7'][np.where(wnoise['pID']==pid)[0]]
        noisegainunscaled=(tsys_all[rca]['tsys_y'][tsys_all[rca]['pID']==pid]+tsky)/wn
        #require calibrated noise to match Tsys*constant (=1/sqrt(bt))
        if len(noisegainunscaled)>0:
            pidout.append(pid)
            noisegainout.append(noisegainunscaled)            
    pidout=np.array(pidout)
    noisegainout=np.array(noisegainout)
    scale=mediannominalgain/nanmedian(noisegainout)
    noisegainout=noisegainout*scale
    noisegainout=removeoutliers(noisegainout.flatten())
    smoothednoisegainout=np.convolve(noisegainout, w, mode='same')
    smoothednoisegainout[:binsize/2]=smoothednoisegainout[binsize/2+1]
    smoothednoisegainout[-binsize/2:]=smoothednoisegainout[-binsize/2-1]
    #make a plot and save it
    plt.hold(False)
    plt.plot(pidout,noisegainout,label='Noise gain')
    plt.hold(True)
    plt.plot(pidout,smoothednoisegainout,label='Binned '+str(binsize),color='black')
    plt.plot(ucds[rca]['nominal_gain'],label='UCDS nominal gain',color='red',lw=2)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.xlabel('Pointing  ID'),plt.ylabel('Gain, K/V'),plt.title(rca + ' Noise estimated gain using Y factor Tsys')
    plt.show()
    plt.savefig('noisegain_'+rca+'.png')
    
    chcalib.append(smoothednoisegainout)
    choff.append(np.zeros(len(smoothednoisegainout)))
write_calfile('tswn',44,pid44,rcalist[13:19],chcalib,choff)


#30 GHz last
freq=30
tsky=.76*2.735
cols=[]
chcalib=[]
choff=[]
#need to find universal intersection PIDs among 70 Ghz radiometers:
pid30=ucds['LFI27M']['pID']
for rca in rcalist[18:]:
    pid30=np.intersect1d(pid30,ucds[rca]['pID'])

for rca in rcalist[18:]:
    #get average nominal gain to match
    mediannominalgain=nanmedian(ucds[rca]['nominal_gain'])
    #deparse to the uncal noise file name
    fnoise=open('wn_'+rca+'.pkl','rb')
    wnoise=cPickle.load(fnoise)
    fnoise.close()
    pidout=[]
    noisegainout=[]
    for pid in pid30:
        wn=wnoise['wn7'][np.where(wnoise['pID']==pid)[0]]
        noisegainunscaled=(tsys_all[rca]['tsys_y'][tsys_all[rca]['pID']==pid]+tsky)/wn
        #require calibrated noise to match Tsys*constant (=1/sqrt(bt))
        if len(noisegainunscaled)>0:
            pidout.append(pid)
            noisegainout.append(noisegainunscaled)            
    pidout=np.array(pidout)
    noisegainout=np.array(noisegainout)
    scale=mediannominalgain/nanmedian(noisegainout)
    noisegainout=noisegainout*scale
    noisegainout=removeoutliers(noisegainout.flatten())
    smoothednoisegainout=np.convolve(noisegainout, w, mode='same')
    smoothednoisegainout[:binsize/2]=smoothednoisegainout[binsize/2+1]
    smoothednoisegainout[-binsize/2:]=smoothednoisegainout[-binsize/2-1]
    #make a plot and save it
    plt.hold(False)
    plt.plot(pidout,noisegainout,label='Noise gain')
    plt.hold(True)
    plt.plot(pidout,smoothednoisegainout,label='Binned '+str(binsize),color='black')
    plt.plot(ucds[rca]['nominal_gain'],label='UCDS nominal gain',color='red',lw=2)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.xlabel('Pointing  ID'),plt.ylabel('Gain, K/V'),plt.title(rca + ' Noise estimated gain using Y factor Tsys')
    plt.show()
    plt.savefig('noisegain_'+rca+'.png')
    
    chcalib.append(smoothednoisegainout)
    choff.append(np.zeros(len(smoothednoisegainout)))
write_calfile('tswn',30,pid30,rcalist[19:],chcalib,choff)


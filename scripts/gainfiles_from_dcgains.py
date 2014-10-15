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
        
#section to include diode weights to compare diode 0 and 1 with appropriate 'nominal gain'

weight_dict={
"LFI18M00" :	0.567304963, 
"LFI18M01" :	0.432695037,
"LFI18S10" :	0.387168785,
"LFI18S11" :	0.612831215,
"LFI19M00" :	0.502457723,
"LFI19M01" :	0.497542277,
"LFI19S10" :	0.55143474,
"LFI19S11" :	0.44856526,
"LFI20M00" :	0.523020094,
"LFI20M01" :	0.476979906,
"LFI20S10" :	0.476730576,
"LFI20S11" :	0.523269424,
"LFI21M00" :	0.500324722,
"LFI21M01" :	0.499675278,
"LFI21S10" :	0.563712153,
"LFI21S11" :	0.436287847,
"LFI22M00" :	0.536283158,
"LFI22M01" :	0.463716842,
"LFI22S10" :	0.553913461,
"LFI22S11" :	0.446086539,
"LFI23M00" :	0.508036034,
"LFI23M01" :	0.491963966,
"LFI23S10" :	0.36160661,
"LFI23S11" :	0.63839339,
"LFI24M00" :	0.602269189,
"LFI24M01" :	0.397730811,
"LFI24S10" :	0.456037835,
"LFI24S11" :	0.543962165,
"LFI25M00" :	0.482050606,
"LFI25M01" :	0.517949394,
"LFI25S10" :	0.369618239,
"LFI25S11" :	0.630381761,
"LFI26M00" :	0.593126369,
"LFI26M01" :	0.406873631,
"LFI26S10" :	0.424268188,
"LFI26S11" :	0.575731812,
"LFI27M00" :	0.519877701,
"LFI27M01" :	0.480122299,
"LFI27S10" :	0.484831449,
"LFI27S11" :	0.515168551,
"LFI28M00" :	0.553227696,
"LFI28M01" :	0.446772304,
"LFI28S10" :	0.467677355,
"LFI28S11" :	0.532322645
}

rj2planck4k={
"LFI18M00" :	0.51, 
"LFI18M01" :	0.68,
"LFI18S10" :	0.68,
"LFI18S11" :	0.68,
"LFI19M00" :	0.68,
"LFI19M01" :	0.68,
"LFI19S10" :	0.68,
"LFI19S11" :	0.68,
"LFI20M00" :	0.68,
"LFI20M01" :	0.68,
"LFI20S10" :	0.68,
"LFI20S11" :	0.68,
"LFI21M00" :	0.68,
"LFI21M01" :	0.68,
"LFI21S10" :	0.68,
"LFI21S11" :	0.68,
"LFI22M00" :	0.68,
"LFI22M01" :	0.68,
"LFI22S10" :	0.68,
"LFI22S11" :	0.68,
"LFI23M00" :	0.68,
"LFI23M01" :	0.68,
"LFI23S10" :	0.68,
"LFI23S11" :	0.68,
"LFI24M00" :	0.79,
"LFI24M01" :	0.79,
"LFI24S10" :	0.79,
"LFI24S11" :	0.79,
"LFI25M00" :	0.79,
"LFI25M01" :	0.79,
"LFI25S10" :	0.79,
"LFI25S11" :	0.79,
"LFI26M00" :	0.79,
"LFI26M01" :	0.79,
"LFI26S10" :	0.79,
"LFI26S11" :	0.79,
"LFI27M00" :	0.85,
"LFI27M01" :	0.85,
"LFI27S10" :	0.85,
"LFI27S11" :	0.85,
"LFI28M00" :	0.85,
"LFI28M01" :	0.85,
"LFI28S10" :	0.85,
"LFI28S11" :	0.85
}

rj2plancksky={
"LFI18M00" :	0.51, 
"LFI18M01" :	0.51,
"LFI18S10" :	0.51,
"LFI18S11" :	0.51,
"LFI19M00" :	0.51,
"LFI19M01" :	0.51,
"LFI19S10" :	0.51,
"LFI19S11" :	0.51,
"LFI20M00" :	0.51,
"LFI20M01" :	0.51,
"LFI20S10" :	0.51,
"LFI20S11" :	0.51,
"LFI21M00" :	0.51,
"LFI21M01" :	0.51,
"LFI21S10" :	0.51,
"LFI21S11" :	0.51,
"LFI22M00" :	0.51,
"LFI22M01" :	0.51,
"LFI22S10" :	0.51,
"LFI22S11" :	0.51,
"LFI23M00" :	0.51,
"LFI23M01" :	0.51,
"LFI23S10" :	0.51,
"LFI23S11" :	0.51,
"LFI24M00" :	0.66,
"LFI24M01" :	0.66,
"LFI24S10" :	0.66,
"LFI24S11" :	0.66,
"LFI25M00" :	0.66,
"LFI25M01" :	0.66,
"LFI25S10" :	0.66,
"LFI25S11" :	0.66,
"LFI26M00" :	0.66,
"LFI26M01" :	0.66,
"LFI26S10" :	0.66,
"LFI26S11" :	0.66,
"LFI27M00" :	0.76,
"LFI27M01" :	0.76,
"LFI27S10" :	0.76,
"LFI27S11" :	0.76,
"LFI28M00" :	0.76,
"LFI28M01" :	0.76,
"LFI28S10" :	0.76,
"LFI28S11" :	0.76
}

#From villa, A&A prelaunch calibration paper
rca_tsys={"LFI18M00" :   36.0 ,	
"LFI18M01" :36.1,	
"LFI18S10" :33.9,	
"LFI18S11" :35.1,
"LFI19M00" :33.1,
"LFI19M01" :31.5,	
"LFI19S10" :32.2,	
"LFI19S11" :33.6,
"LFI20M00" :35.2,
"LFI20M01" :34.2,	
"LFI20S10" :36.9,	
"LFI20S11" :35.0,
"LFI21M00" :27.3,
"LFI21M01" :28.4,	
"LFI21S10" :34.4,	
"LFI21S11" :36.4,
"LFI22M00" :30.9,
"LFI22M01" :30.3,	
"LFI22S10" :30.3,	
"LFI22S11" :31.8,
"LFI23M00" :35.9,
"LFI23M01" :34.1,	
"LFI23S10" :33.9,	
"LFI23S11" :31.1,
"LFI24M00" :15.5,
"LFI24M01" :15.3,	
"LFI24S10" :15.8,	
"LFI24S11" :15.8,
"LFI25M00" :17.5,
"LFI25M01" :17.9,	
"LFI25S10" :18.6,	
"LFI25S11" :18.4,
"LFI26M00" :18.4,
"LFI26M01" :17.4,	
"LFI26S10" :16.8,	
"LFI26S11" :16.5,
"LFI27M00" :12.1,
"LFI27M01" :11.9,	
"LFI27S10" :13.0,	
"LFI27S11" :12.5,
"LFI28M00" :10.6,
"LFI28M01" :10.3,	
"LFI28S10" :9.9,
"LFI28S11" :9.8}

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

ucds30=pu.get_ucds(freq=30)
ucds44=pu.get_ucds(freq=44)
ucds70=pu.get_ucds(freq=70)
ucds=dict(ucds70,**ucds44)
ucds.update(ucds30)
#set up the output moving average:
binsize=30
w = np.empty((binsize,), dtype=np.float_)
w[:] = 1.0/binsize



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
for rca in rcalist[:12]:
    pid70=np.intersect1d(pid70,ucds[rca]['pID'])

for key in rcalist[:12]:
    if key[-1]=='S':
        diode='10'
        tsys_rca=(rca_tsys[key+'10'] + rca_tsys[key+'11'])/2.0
    if key[-1]=='M':
        diode='00'
        tsys_rca=(rca_tsys[key+'00'] + rca_tsys[key+'01'])/2.0
    t4k=ucds[key]['l14k']
    t4k[t4k==0]=4.655
    thot=rj2planck4k[key+diode]*t4k
    tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(ucds[key])))
    offr=nanmedian(ucds[key]['vref']-(tsys_rca+thot)/ucds[key]['nominal_gain'])
    offs=nanmedian(ucds[key]['vsky']-(tsys_rca+tcold)/ucds[key]['nominal_gain'])
    tsys_ref=(ucds[key]['vref']-offr)*ucds[key]['nominal_gain']-thot
    tsys_sky=(ucds[key]['vsky']-offs)*ucds[key]['nominal_gain']-tcold
    dcgain=(thot- tcold)/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
    dcgain=removeoutliers(dcgain)
    scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(dcgain)
    dcgain=dcgain*scale
    dcgainout=[]
    for pid in pid70:
        dcgainout.append(dcgain[ucds[key]['pID']==pid])
    pidout=np.array(pid70).flatten()
    dcgainout=np.array(dcgainout).flatten()
    
    smootheddcgainout=np.convolve(dcgainout, w, mode='same')
    smootheddcgainout[:binsize/2]=smootheddcgainout[binsize/2+1]
    smootheddcgainout[-binsize/2:]=smootheddcgainout[-binsize/2-1]
    #make a plot and save it
    plt.hold(False)
    plt.plot(pidout,dcgainout,label='DC gain R-S')
    plt.hold(True)
    plt.plot(pidout,smootheddcgainout,label='Binned '+str(binsize),color='black')
    plt.plot(ucds[rca]['nominal_gain'],label='UCDS nominal gain',color='red',lw=2)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.xlabel('Pointing  ID'),plt.ylabel('Gain, K/V'),plt.title(rca + 'UCDS dc estimated gain using Ref-Sky')
    plt.show()
    plt.savefig('dcgain_'+rca+'.png')
    
    chcalib.append(dcgainout)
    choff.append(np.zeros(len(dcgainout)))
write_calfile('DCRS',70,pid70,rcalist[:12],chcalib,choff)


#44 GHz next
freq=44
tsky=.66*2.735
cols=[]
chcalib=[]
choff=[]

#need to find universal intersection PIDs among 70 Ghz radiometers:
pid44=ucds['LFI24M']['pID']
for rca in rcalist[12:18]:
    pid44=np.intersect1d(pid44,ucds[rca]['pID'])

for key in rcalist[12:18]:
    if key[-1]=='S':
        diode='10'
        tsys_rca=(rca_tsys[key+'10'] + rca_tsys[key+'11'])/2.0
    if key[-1]=='M':
        diode='00'
        tsys_rca=(rca_tsys[key+'00'] + rca_tsys[key+'01'])/2.0
    t4k=ucds[key]['l14k']
    t4k[t4k==0]=4.655
    thot=rj2planck4k[key+diode]*t4k
    tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(ucds[key])))
    offr=nanmedian(ucds[key]['vref']-(tsys_rca+thot)/ucds[key]['nominal_gain'])
    offs=nanmedian(ucds[key]['vsky']-(tsys_rca+tcold)/ucds[key]['nominal_gain'])
    tsys_ref=(ucds[key]['vref']-offr)*ucds[key]['nominal_gain']-thot
    tsys_sky=(ucds[key]['vsky']-offs)*ucds[key]['nominal_gain']-tcold
    dcgain=(thot- tcold)/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
    dcgain=removeoutliers(dcgain)
    scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(dcgain)
    dcgain=dcgain*scale
    dcgainout=[]
    for pid in pid44:
        dcgainout.append(dcgain[ucds[key]['pID']==pid])
    pidout=np.array(pid44).flatten()
    dcgainout=np.array(dcgainout).flatten()
    
    smootheddcgainout=np.convolve(dcgainout, w, mode='same')
    smootheddcgainout[:binsize/2]=smootheddcgainout[binsize/2+1]
    smootheddcgainout[-binsize/2:]=smootheddcgainout[-binsize/2-1]
    #make a plot and save it
    plt.hold(False)
    plt.plot(pidout,dcgainout,label='DC gain R-S')
    plt.hold(True)
    plt.plot(pidout,smootheddcgainout,label='Binned '+str(binsize),color='black')
    plt.plot(ucds[rca]['nominal_gain'],label='UCDS nominal gain',color='red',lw=2)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.xlabel('Pointing  ID'),plt.ylabel('Gain, K/V'),plt.title(rca + 'UCDS dc estimated gain using Ref-Sky')
    plt.show()
    plt.savefig('dcgain_'+rca+'.png')
    
    chcalib.append(dcgainout)
    choff.append(np.zeros(len(dcgainout)))
write_calfile('DCRS',44,pid44,rcalist[12:18],chcalib,choff)


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

for key in rcalist[18:]:
    if key[-1]=='S':
        diode='10'
        tsys_rca=(rca_tsys[key+'10'] + rca_tsys[key+'11'])/2.0
    if key[-1]=='M':
        diode='00'
        tsys_rca=(rca_tsys[key+'00'] + rca_tsys[key+'01'])/2.0
    t4k=ucds[key]['l14k']
    t4k[t4k==0]=4.655
    thot=rj2planck4k[key+diode]*t4k
    tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(ucds[key])))
    offr=nanmedian(ucds[key]['vref']-(tsys_rca+thot)/ucds[key]['nominal_gain'])
    offs=nanmedian(ucds[key]['vsky']-(tsys_rca+tcold)/ucds[key]['nominal_gain'])
    tsys_ref=(ucds[key]['vref']-offr)*ucds[key]['nominal_gain']-thot
    tsys_sky=(ucds[key]['vsky']-offs)*ucds[key]['nominal_gain']-tcold
    dcgain=(thot- tcold)/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
    dcgain=removeoutliers(dcgain)
    scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(dcgain)
    dcgain=dcgain*scale
    dcgainout=[]
    for pid in pid30:
        dcgainout.append(dcgain[ucds[key]['pID']==pid])
    pidout=np.array(pid30).flatten()
    dcgainout=np.array(dcgainout).flatten()
    
    smootheddcgainout=np.convolve(dcgainout, w, mode='same')
    smootheddcgainout[:binsize/2]=smootheddcgainout[binsize/2+1]
    smootheddcgainout[-binsize/2:]=smootheddcgainout[-binsize/2-1]
    #make a plot and save it
    plt.hold(False)
    plt.plot(pidout,dcgainout,label='DC gain R-S')
    plt.hold(True)
    plt.plot(pidout,smootheddcgainout,label='Binned '+str(binsize),color='black')
    plt.plot(ucds[rca]['nominal_gain'],label='UCDS nominal gain',color='red',lw=2)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.xlabel('Pointing  ID'),plt.ylabel('Gain, K/V'),plt.title(rca + 'UCDS dc estimated gain using Ref-Sky')
    plt.show()
    plt.savefig('dcgain_'+rca+'.png')
    
    chcalib.append(dcgainout)
    choff.append(np.zeros(len(dcgainout)))
write_calfile('DCRS',30,pid30,rcalist[18:],chcalib,choff)




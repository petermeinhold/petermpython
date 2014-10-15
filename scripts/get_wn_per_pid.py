"""
script to read in DPC data and estimate uncalibrated white noise for each PID
"""
from glob import glob
from scipy.stats import nanmean,nanstd,nanmedian
import pyfits
import cPickle
import planck_util as pu

#first read in all the ucds so we can get the PID OBT starts and stops
ucds30=pu.get_ucds(freq=30)
ucds44=pu.get_ucds(freq=44)
ucds70=pu.get_ucds(freq=70)
ucds=dict(ucds70,**ucds44)
ucds.update(ucds30)

def get_wn_per_pid(rca='LFI18M'):
    horn=int(rca[3:5])
    freq=30
    samprate=32.51
    if horn<27:
        freq=44
        samprate=46.55
    if horn<24:
        freq=70
        samprate=78.76
    # start with pidlist,obtstart,obtend and odlist
    pidlist=ucds[rca]['pID']
    obtstarts=ucds[rca]['obt_start']
    obtends=ucds[rca]['obt_end']
    ods=ucds[rca]['od_int']
    odlist=unique(ods)
    odlist.sort()
    wnlist9=[]
    wnlist7=[]
    wnlist5=[]
    pidoutlist=[]
    obtoutlist=[]
    for od in odlist:
        print od
        data=pu.get_lfi_dpc_timeline(od,rca)
        for pid in pidlist[ods==od]:
            piddata=data[1][(data[0] < obtends[pidlist==pid][0]) & (data[0] > obtstarts[pidlist==pid][0])]
            if len(piddata)>20.*samprate:
                z=pu.nps(piddata,samprate,minfreq=.1)
                wnlist9.append(sqrt(np.mean(z[1][z[0] > .9*samprate/2.])))
                wnlist7.append(sqrt(np.mean(z[1][z[0] > .7*samprate/2.])))
                wnlist5.append(sqrt(np.mean(z[1][z[0] > .5*samprate/2.])))
                pidoutlist.append(pid)
                obtoutlist.append(obtstarts[pidlist==pid][0])
    return {'wn9':wnlist9,'wn7':wnlist7,'wn5':wnlist5,'pID':pidoutlist,'obtstart':obtoutlist}

wn_channels18_19_20={}
for r in range(18,21):
    for d in ['M','S']:
        rca='LFI'+str(r)+d
        print rca
        wn=get_wn_per_pid(rca=rca)
        wn_channels18_19_20.update({rca:wn})
f=open('wn_uncal_channels18_19_20.pkl','wb')
cPickle.dump(wn_channels18_19_20,f)
f.close()

wn_channels21_22_23={}
for r in range(21,24):
    for d in ['M','S']:
        rca='LFI'+str(r)+d
        print rca
        wn=get_wn_per_pid(rca=rca)
        wn_channels21_22_23.update({rca:wn})
f=open('wn_uncal_channels21_22_23.pkl','wb')
cPickle.dump(wn_channels21_22_23,f)
f.close()

wn_channels24_25_26_27_28={}
for r in range(24,29):
    for d in ['M','S']:
        rca='LFI'+str(r)+d
        print rca
        wn=get_wn_per_pid(rca=rca)
        wn_channels24_25_26_27_28.update({rca:wn})
f=open('wn_uncal_channels24_25_26_27_28.pkl','wb')
cPickle.dump(wn_channels24_25_26_27_28,f)
f.close()




            
            
            
""" 
multiprocess test script, try running all  70 GHz with 8 processes
script to read in DPC data and estimate uncalibrated white noise for each PID
"""
from glob import glob
from scipy.stats import nanmean,nanstd,nanmedian
import pyfits
import cPickle
import planck_util as pu
from multiprocessing import Pool
import numpy as np

#first read in all the ucds so we can get the PID OBT starts and stops


def get_wn_per_pid(rca):
    ucds30=pu.get_ucds(freq=30)
    ucds44=pu.get_ucds(freq=44)
    ucds70=pu.get_ucds(freq=70)
    ucds=dict(ucds70,**ucds44)
    ucds.update(ucds30)
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
    odlist=np.unique(ods)
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
                wnlist9.append(np.sqrt(np.mean(z[1][z[0] > .9*samprate/2.])))
                wnlist7.append(np.sqrt(np.mean(z[1][z[0] > .7*samprate/2.])))
                wnlist5.append(np.sqrt(np.mean(z[1][z[0] > .5*samprate/2.])))
                pidoutlist.append(pid)
                obtoutlist.append(obtstarts[pidlist==pid][0])
    wnoutput= {'wn9':wnlist9,'wn7':wnlist7,'wn5':wnlist5,'pID':pidoutlist,'obtstart':obtoutlist}
    f=open('wn_'+rca+'.pkl','wb')
    cPickle.dump(wnoutput,f)
    f.close()


if __name__=='__main__':

    rcalist=[]
    for r in range(18,24):
        for d in ['M','S']:
            rcalist.append('LFI'+str(r)+d)
    
    pool=Pool(processes=12)
    pool.map(get_wn_per_pid,rcalist)
    





            
            
            
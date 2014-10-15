import planck_util as pu
from glob import glob
import pyfits
from testenv.ringset import *
import healpy as hp
from scipy.stats import nanmean,nanstd,nanmedian
from testenv.todtools import write_calfile

#load the  ucds just to make OD to pID lookup table:

ucds=pu.get_ucds(freq=30)


def od2pid_range(od,ucds=None,key='LFI27M'):
    """
    function returns 2 element range min and max pID of the input OD. If UCDS is None, reads it in
    """
    if ucds==None:
        ucds=pu.get_ucds(freq=30)
    pidrange=[np.min(ucds[key]['pID'][ucds[key]['od_int']==od]),np.max(ucds[key]['pID'][ucds[key]['od_int']==od])]
    return pidrange

def surv_diff(surveydict,ss1=1,ss2=2,nside=128,freq=70,mask_ps=True,fwhm=10.0):
    """
    function to make differences among 5 surveys from andreas interactive destriper/binner. uses common masks smooth to 10 degrees
    fwhm is in degrees, if zero or negative don't smooth at all, default 10 degrees
    """
    m1=hp.ma(np.array(surveydict[ss1]))
    m2=hp.ma(np.array(surveydict[ss2])) 
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%GHz_*.fits' % freq)[0]), nside,order_out='NEST')))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    dif.mask |= np.isnan(dif)
    dif=hp.ud_grade(dif,nside,order_in='NEST',order_out='RING')
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=dif.mask
    if fwhm<=0:
        difsm=dif
    return difsm
    
def surv_sum_diff(surveydict,ss=1,nside=128,freq=70,mask_ps=True,fwhm=10.0):
    """
    function to genreate smoothed difference between chosen survey and sum of the other 5 surveys from andreas interactive destriper/binner. uses common masks smooth to 10 degrees
    fwhm is in degrees, if zero or negative don't smooth at all, default 10 degrees
    """
    sumlist=surveydict.keys()
    sumlist.remove(ss)
    nsum=len(sumlist)
    m=hp.ma(np.array(surveydict[sumlist[0]]))/nsum
    totalmask=m.mask
    if mask_ps==True:
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside,order_out='NEST')))
        totalmask=m.mask | psmask
    for ss in sumlist[1:]:
        m1=hp.ma(np.array(surveydict[ss]))
        m=m+m1/nsum
        totalmask=m1.mask | totalmask
    m.mask=totalmask
    m.mask |= np.isnan(m)
    
    m1=hp.ma(np.array(surveydict[ss]))
    totalmask=totalmask | m1.mask
    d=(m1-m)/2.
    d=hp.ud_grade(d,nside,order_in='NEST',order_out='RING')
    if fwhm>0:
        dsm=hp.ma(hp.smoothing(d.filled(),fwhm*np.pi/180.))
        dsm.mask=m.mask
    if fwhm<=0:
        dsm=d
    hp.mollview(dsm,min=-1e-5,max=1e-5,title='SS'+np.str(ss)+ '  - sum of others')
    return dsm

def bin_test_cal_to_od(caltype='tswn',freq=70,rca='LFI19M',relative=True,odlist=None):
    #function to read in candidate calibration and bin to od to use in mapmaker
    print caltype
    if caltype != 'nominal':
        caldir='/project/projectdirs/planck/user/zonca/cal/'+caltype+'/'
        f=glob(caldir+'*'+str(freq)+'*.fits')[-1]
        ff=pyfits.open(f)
        gain_pid=ff[rca].data[rca]
        pid=ff['pid'].data['PID']
    ucds=pu.get_ucds(freq=freq)
    if caltype=='nominal':
        gain_pid=ucds[rca]['nominal_gain']
        pid=ucds[rca]['pID']
    if odlist==None:
        odlist=np.unique(ucds[rca]['od_int'])
    odlist.sort()
    
    gain=[]
    odoutlist=[]
    gain_nom=[]
    for od in odlist:
        u_odpids=ucds[rca]['pID'][ucds[rca]['od_int']==od]
        if len(u_odpids)>0:   #if we dont have any PID's for that OD, use the last ones
            minpid=np.min(u_odpids)
            maxpid=np.max(u_odpids)
            odgain_nom=np.mean(ucds[rca]['nominal_gain'][ucds[rca]['od_int']==od])
        pidod=where((pid>=minpid) & (pid <= maxpid))
        odgain=np.mean(gain_pid[pidod])

        if odgain > 0:
            gain.append(odgain)
            odoutlist.append(od)
            gain_nom.append(odgain_nom)
    
#    if caltype=='nominal':
#        gain=np.ones(len(odlist))
#        odoutlist=odlist
#        relative=False
#        print 'in nominal, set to 1'
    gain=np.array(gain)
    odoutlist=np.array(odoutlist)
    gain_nom=np.array(gain_nom)
    plt.figure()
    plt.plot(odoutlist,gain,label=caltype+' gain')
    plt.xlabel('OD'),plt.ylabel('Calibration, K/V'),plt.title(caltype + '  Calibration per OD '+rca)
    plt.plot(odoutlist,gain_nom,label='nominal gain')
    plt.ylim([np.mean(gain)-4*np.std(gain),np.mean(gain)+4*np.std(gain)])
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.savefig('/global/homes/p/peterm/plots/cal_comparisons/cal_'+caltype+'_'+rca+'.png')
    if relative:
        gain=gain/np.mean(gain_nom)
    return odoutlist,gain
    
    
    
def caltest_surveys(caltype='live',freq=70,rca='LFI19M',nside=128,pids=None,cals=None,calname=None):
    R = RingSetManager(rca, nside)
    #R.apply_mask("destripingmask_30.fits")
    ods = ods_from_tag("full")
    hit_map = R.create_hit_map(ods)
    if caltype != 'live':
        odcal,calibration = bin_test_cal_to_od(caltype=caltype,odlist=ods,rca=rca)
    if caltype == 'live':
        odcal,calibration=bin_test_cal_to_od_live(pids,cals,rca=rca,freq=freq,calname=calname)
    plt.plot(odcal,calibration)

    cal = pd.Series(calibration, index=odcal)

    calibrated_ringsets = R.calibrate(cal)
    bin_map, destriped_map, baseline_removed = R.destripe(calibrated_ringsets, return_baseline_removed=True)

    survey_map = {}
    for survey in range(1, 5+1):
        print "Binning survey %d" % survey
        survey_map[survey] = R.create_bin_map(slice_data(baseline_removed, "survey%d" % survey))
        #hp.mollview(survey_map[survey], nest=True, min=-1e-3, max=1e-3, title="%s survey %d" % (R.ch.tag, survey))
    return survey_map
    
    
def subtract(map1,map2):
    map1=hp.ma(map1)
    map2=hp.ma(map2)
    diff=hp.ma((map1-map2)/2.)
    diff.mask=map1.mask | map2.mask
    return diff
    
def dif_map(rca='LFI19M',nside=128):
    """
    script to make maps and survey diffs smoothed 10 degrees. save relative gain plot and maps to disk
    """
    freq=70
    h=int(rca[3:5])
    if h>23:
        freq=44
    if h>26:
        freq=30

    caltypes=['nominal']  #,'gain1'] #,'tswn','DCRS']
    for caltype in caltypes:
        s_surveys=caltest_surveys(caltype=caltype,freq=freq,rca=rca)
        for ss1,ss2 in zip([1,1,1,1],[2,3,4,5]):
            d=surv_diff(s_surveys,ss1=ss1,ss2=ss2,nside=128,freq=freq,mask_ps=True,fwhm=10.0)
            hp.mollview(d,min=-10e-6,max=10e-6,title=rca+'  '+caltype+'   '+ 'SS'+str(ss1) + ' - '+'SS'+str(ss2))
            savefig('diffmap_'+rca+'_'+caltype+'_'+ 'SS'+str(ss1) + '-'+'SS'+str(ss2)+'.png')
    


def dif_map_live(rca='LFI19M',nside=128,pids=None,cal_list=None,cal_names=None):
    """
    script to make maps and survey diffs smoothed 10 degrees. save relative gain plot and maps to disk
    """
    freq=70
    h=int(rca[3:5])
    if h>23:
        freq=44
    if h>26:
        freq=30

    for calname,cals in zip(cal_names,cal_list):
        s_surveys=caltest_surveys(caltype='live',freq=freq,rca=rca,pids=pids,cals=cals,calname=calname)
        for ss1,ss2 in zip([1,1,1,1],[2,3,4,5]):
            d=surv_diff(s_surveys,ss1=ss1,ss2=ss2,nside=128,freq=freq,mask_ps=True,fwhm=10.0)
            hp.mollview(d,min=-10e-6,max=10e-6,title=rca+'  '+calname +'   '+ 'SS'+str(ss1) + ' - '+'SS'+str(ss2))
            savefig('diffmap_'+rca+'_'+calname+'_'+ 'SS'+str(ss1) + '-'+'SS'+str(ss2)+'.png')
    


def bin_test_cal_to_od_live(inputpids,cals,freq=70,rca='LFI19M',relative=True,odlist=None,calname='live'):
    #function to bin 'live' test calibration to odlist
    gain_pid=cals
    pid=inputpids
    ucds=pu.get_ucds(freq=freq)
    if odlist==None:
        odlist=np.unique(ucds[rca]['od_int'])
    odlist.sort()
    
    gain=[]
    odoutlist=[]
    gain_nom=[]
    for od in odlist:
        u_odpids=ucds[rca]['pID'][ucds[rca]['od_int']==od]
        if len(u_odpids)>0:   #if we dont have any PID's for that OD, use the last ones
            minpid=np.min(u_odpids)
            maxpid=np.max(u_odpids)
            odgain_nom=np.mean(ucds[rca]['nominal_gain'][ucds[rca]['od_int']==od])
        pidod=where((pid>=minpid) & (pid <= maxpid))
        odgain=np.mean(gain_pid[pidod])

        if odgain > 0:
            gain.append(odgain)
            odoutlist.append(od)
            gain_nom.append(odgain_nom)
    
#    if caltype=='nominal':
#        gain=np.ones(len(odlist))
#        odoutlist=odlist
#        relative=False
#        print 'in nominal, set to 1'
    gain=np.array(gain)
    odoutlist=np.array(odoutlist)
    gain_nom=np.array(gain_nom)
    plt.figure()
    plt.plot(odoutlist,gain,label=calname+' gain')
    plt.xlabel('OD'),plt.ylabel('Calibration, K/V'),plt.title(calname + '  Calibration per OD '+rca)
    plt.plot(odoutlist,gain_nom,label='nominal gain')
    plt.ylim([np.mean(gain)-4*np.std(gain),np.mean(gain)+4*np.std(gain)])
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.savefig('/global/homes/p/peterm/plots/cal_comparisons/cal_'+calname+'_'+rca+'.png')
    if relative:
        gain=gain/np.mean(gain_nom)
    return odoutlist,gain

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



ucds30=pu.get_ucds(freq=30)
ucds44=pu.get_ucds(freq=44)
ucds70=pu.get_ucds(freq=70)
ucds=dict(ucds70,**ucds44)
ucds.update(ucds30)


rca='LFI19M'
key=rca


if key[-1]=='S':
    diode='10'
    tsys_rca=(rca_tsys[key+'10'] + rca_tsys[key+'11'])/2.0
if key[-1]=='M':
    diode='00'
    tsys_rca=(rca_tsys[key+'00'] + rca_tsys[key+'01'])/2.0
t4k=ucds[key]['l14k']
t4k[t4k==0]=4.655

t4k=ucds[key]['l14k']
t4k[t4k==0]=4.655
thot=rj2planck4k[key+diode]*t4k
tfem=ucds[key]['lfem1']
tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(ucds[key])))
offr=nanmedian(ucds[key]['vref']-(tsys_rca+thot)/ucds[key]['nominal_gain'])
offs=nanmedian(ucds[key]['vsky']-(tsys_rca+tcold)/ucds[key]['nominal_gain'])


rcalist=[]
for r in range(18,29):
    for d in ['M','S']:
        rcalist.append('LFI'+str(r)+d)
        

#gain_nominal:
gain=(ucds[key]['nominal_gain'])
gain=removeoutliers(gain)
gain_nominal=gain[1:]
pids=ucds[rca]['pID'][1:]

   
#gain1:
lc=.985
lh=.995 
offr=0.19 
offs=0.2
gain=(lh*thot+(1-lh)*tfem- (lc*tcold+(1-lc)*tfem))/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
gain=removeoutliers(gain)
scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(gain)
print scale
gain1=gain[1:]*scale
pids=ucds[rca]['pID'][1:]

#gain2:
lc=.985
lh=.995 
offr=0.19 
offs=0.3
gain=(lh*thot+(1-lh)*tfem- (lc*tcold+(1-lc)*tfem))/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
gain=removeoutliers(gain)
scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(gain)
print scale
gain2=gain[1:]*scale
pids=ucds[rca]['pID'][1:]

#gain3:
lc=.89
lh=.89
offr=0.1 
offs=0.19
gain=(lh*thot+(1-lh)*tfem- (lc*tcold+(1-lc)*tfem))/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
gain=removeoutliers(gain)
scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(gain)
print scale
gain3=gain[1:]*scale
pids=ucds[rca]['pID'][1:]









for a in adj:
    offr=offr+a
    lh=.95
    lc=.95
    gain=(lh*thot+(1-lh)*tfem- (lc*tcold+(1-lc)*tfem))/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
    gain=removeoutliers(gain)
    scale=nanmedian(ucds[key]['nominal_gain'])/nanmedian(gain)
    print scale
    gain=gain[1:]*scale
    pids=ucds[rca][pID][1:]
    plot ucds[rca]['pID'][1:],gain[1:],label='lc:'+str(lc)+'  lh:'+str(lh) +'  or:'+str(offr)[:4]+'  os:'+str(offs)[:4]
    
    
    trygains.append(gain)
    
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
    
  

ladj=arange(.8,1.0,.02)
oadj=arange(.1,.2,.01)
sigs=[]
offss=[]
offrs=[]
lcs=[]
lhs=[]
nom_median=nanmedian(ucds[key]['nominal_gain'])
for lc in ladj:
    for lh in ladj:
        for offr in oadj:
            for offs in oadj:
                gain=(lh*thot+(1-lh)*tfem- (lc*tcold+(1-lc)*tfem))/((ucds[key]['vref']-offr)-(ucds[key]['vsky']-offs))
                gain=removeoutliers(gain)
                scale=nom_median/nanmedian(gain)
                gain=gain*scale
                sigs.append(np.std(gain-ucds[key]['nominal_gain']))
                lcs.append(lc)
                lhs.append(lh)
                offss.append(offs)
                offrs.append(offr)
                



    
    
def example():
    R = RingSetManager("LFI19M", 128)
    #R.apply_mask("destripingmask_30.fits")
    ods = ods_from_tag("full")
    hit_map = R.create_hit_map(ods)
    calibration = np.ones(len(ods), dtype=np.float)
    cal = pd.Series(calibration, index=ods)

    calibrated_ringsets = R.calibrate(cal)
    bin_map, destriped_map, baseline_removed = R.destripe(calibrated_ringsets, return_baseline_removed=True)

    survey_map = {}
    
    for survey in range(1, 5+1):
        print "Binning survey %d" % survey
        survey_map[survey] = R.create_bin_map(slice_data(baseline_removed, "survey%d" % survey))
    return survey_map
    #    hp.mollview(survey_map[survey], nest=True, min=-1e-3, max=1e-3, title="%s survey %d" % (R.ch.tag, survey))
    
    
    

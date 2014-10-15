import quickring as qr
from glob import glob
import healpy as hp
import get_rings_pix
import cPickle
from petermpython import mpfit as mp
from petermpython import planck_util as pu

hfiringstart=237
nside=1024
freq=70
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
tmask=wmask.astype(np.bool)|psmask
ncals=5 #number of calibration adjustment periods to use per survey
surveystarts=[3,5484,10958,16455,21483]
surveyends=[5483,10957,16454,21482,27404]


def gaussian_resid(p, fjac=None, x=None, y=None, err=None):
    # Parameter values are passed in "p"
    # for gaussian p=[amplitude,sigma]
    # form is f(x)=p[0]*exp(-((x-p[1])/p[2])^2) 
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.su
    # p=np.array(p)
    model = p[0]*np.exp(-((x-p[1])/p[2])**2)    
    
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return([status, (y-model)/err])
    
def gaussian_model(p,x=None):
    # Parameter values are passed in "p"
    # for gaussian p=[amplitude,sigma]
    # form is f(x)=p[0]*exp(-((x-p[1])/p[2])^2) 
    
    model = p[0]*np.exp(-((x-p[1])/p[2])**2)    
    
    return(model)
    
def fit_a_peak(cl,maxiter=10,yerr=None,frac=None,lmin=100,lmax=350):
    # procedure to run mpfit on a cl, coorect by frac, scale buy l(l+1) and
    # fit from lmin to lmax, assume Cl input is l by l)
    l=np.arange(len(cl))
    if frac==None:
        frac=1.0
    y=l*(l+1)*cl*1e12/(2*np.pi*frac)
    if yerr==None:
        yerr=np.std(y[400:600])
    err=yerr*(np.ones(y.size,dtype=float))
    fa={'x':l[lmin:lmax],'y':y[lmin:lmax],'err':err[lmin:lmax]}
    
    pstart=np.array([4500.,200.,100.])
    m=mp.mpfit(gaussian_resid,pstart,functkw=fa,quiet=1,maxiter=10)
    redchi=m.fnorm
    if m.dof>0:
        redchi=redchi/m.dof
    return(m,redchi)



# this section produced the subsurvey mask
for surv in range(5):
    calstep=np.int(surveyends[surv]-surveystarts[surv])/ncals
    startring=hfiringstart+surveystarts[surv]
    for n in range(ncals-1):
        pxmap=np.ones(hp.nside2npix(1024))
        px=get_rings_pix.get_rings_pix(startring+n*calstep,startring+(n+1)*calstep,nphi=50000)
        pxmap[px]=0
        ttmask=pxmap.astype(np.bool)|tmask
        f=open('subsurveymasks_5/survey_'+np.str(surv+1)+'_sub_'+np.str(n+1)+'.pkl','wb')
        cPickle.dump(ttmask,f)
        f.close()


# new idea, use channel full mission map for cross spectra
base='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/Single_Radiometer/'
basefreq='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/' 
fmapfull=glob(basefreq+'*11_full_*.fits')
m70=hp.ma(hp.read_map(fmapfull[0]))
xpeaks={}
chans=['18M','18S','19M','19S','20M','20S','21M','21S','22M','22S','23M','23S']
uchans=['LFI18M','LFI18S','LFI19M','LFI19S','LFI20M','LFI20S','LFI21M','LFI21S','LFI22M','LFI22S','LFI23M','LFI23S']

for ch in chans:
    xpeaks[ch]={}
    dopp_peaks=[]
    dopp_peak_sigs=[]
    dopp_peak_chis=[]
    doppmeans=[]
    doppsigs=[]
    calpids=[]
    l=arange(3072)
    for surv in range(1,6):
        calstep=np.int(surveyends[surv-1]-surveystarts[surv-1])/ncals
        fl1=glob(base+'*'+ch+'*_survey_'+np.str(surv)+'*.fits')
        m1=hp.ma(hp.read_map(fl1[0]))
        m1mask=copy(m1.mask)
        for n in range(ncals-1):
            f=open('subsurveymasks/survey_'+np.str(surv)+'_sub_'+np.str(n+1)+'.pkl','rb')
            survmask=cPickle.load(f)
            m1.mask=m1mask|survmask
            m70.mask=m1.mask
            cl=hp.anafast(m1,m70,lmax=1535)
            frac=(1-float(m1.mask.sum())/len(m1.mask))
            calpids.append(hfiringstart+surveystarts[surv-1]+n*calstep)
            m,rchi=fit_a_peak(cl,frac=frac)
            dopp_peaks.append(m.params[0])
            dopp_peak_sigs.append(m.perror[0])
            dopp_peak_chis.append(rchi)
            doppmeans.append(np.mean(l[50:250]*((l[50:250]+1)*cl[50:250])/frac))
            doppsigs.append((np.std(l[50:250]*((l[50:250]+1)*cl[50:250])/frac))/np.sqrt(200.))
        doppmeans
        doppsigs
    xpeaks[ch]['fitpeaks']=dopp_peaks
    xpeaks[ch]['fitpeaksigs']=dopp_peak_sigs
    xpeaks[ch]['fitchis']=dopp_peak_chis
    xpeaks[ch]['calpids']=calpids            
    xpeaks[ch]['peakmeans']=doppmeans
    xpeaks[ch]['peaksigs']=doppsigs
f=open('subsurveymasks_5/xpeak_fit_summary.pkl','wb')
cPickle.dump(xpeaks,f)
f.close()


#for n in range(9):
#    plt.figure()
#    for ch in chans:
#        errorbar(float(ch[0:2]),xpeaks[ch]['fitpeaks'][n],yerr=xpeaks[ch]['fitpeaksigs'][n],fmt=None,label=ch)
#    plt.xlabel('Horn number')
#    plt.ylabel('Doppler Peak fit amplitude, microK^2')
#    plt.title('Survey 1, sub '+np.str(n)+'  Chan '+ch)
#    plt.xlim([17,24])
#    leg=plt.legend()
#    plt.savefig('plots/dopp_comparison'+ch+np.str(n)+'.png')
    
ucds=pu.get_ucds(freq=70)
fl=glob('subsurveymasks_5/x*.pkl')
f=open(fl[0],'rb')
xpeaks=cPickle.load(f)
f.close()


#dopp=[]
#for n in range(45):
#    for ch in chans:
#        dopp.append(xpeaks[ch]['fitpeaks'][n])
#dopp=np.array(dopp)
#dopp=dopp.reshape(45,12)
#doppm=dopp-np.mean(dopp,axis=1).reshape(45,1)
#dopp_pids=xpeaks['18M']['calpids']
for ch in chans:
    xpeaks[ch]['normed_gains']=np.ones(len(xpeaks[ch]['fitpeaks']))
for n in range(len(xpeaks[ch]['fitpeaks'])):
    nmean=np.mean([xpeaks[ch]['fitpeaks'][n] for ch in chans])
    for ch in chans:
        xpeaks[ch]['normed_gains'][n]=xpeaks[ch]['fitpeaks'][n]/nmean


testgg={}
for ch,uch in zip(chans,uchans):
    testgg[ch]={}
    upids=ucds[uch]['pID']
    dpids=xpeaks[ch]['calpids']
    ugains=ucds[uch]['nominal_gain']
    dngains=xpeaks[ch]['normed_gains']
    idngains=np.interp(upids,dpids,dngains)
    testgg[ch]['gain']=ugains*idngains
    testgg[ch]['pID']=upids    
    
for c in range(12):
    plt.figure()
    plt.plot(ucds[uchans[c]]['pID'],ucds[uchans[c]]['nominal_gain'] ,label='DPC gain')
    plt.plot(testgg[chans[c]]['pID'],testgg[chans[c]]['gain'] ,label='Modified gain')
    plt.xlabel('PID'),plt.ylabel('Calibration, (K/V)'),plt.title('Gain modification from 5 subsurvey doppler peak estimates '+uchans[c])
    plt.show()
    plt.savefig('subsurveymasks_5/gain'+uchans[c])
    
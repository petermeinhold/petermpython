import cPickle
import healpy as hp
import planck_util as pu
def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()
l28s12=pu.get_surv_diff_chan(chan1='28S',surv1=1,chan2='28S',surv2=2,nside=128)
l28s34=pu.get_surv_diff_chan(chan1='28S',surv1=3,chan2='28S',surv2=4,nside=128)

l27s12=pu.get_surv_diff_chan(chan1='27S',surv1=1,chan2='27S',surv2=2,nside=128)
l27s34=pu.get_surv_diff_chan(chan1='27S',surv1=3,chan2='27S',surv2=4,nside=128)

l28m12=pu.get_surv_diff_chan(chan1='28M',surv1=1,chan2='28M',surv2=2,nside=128)
l28m34=pu.get_surv_diff_chan(chan1='28M',surv1=3,chan2='28M',surv2=4,nside=128)

l27m12=pu.get_surv_diff_chan(chan1='27M',surv1=1,chan2='27M',surv2=2,nside=128)
l27m34=pu.get_surv_diff_chan(chan1='27M',surv1=3,chan2='27M',surv2=4,nside=128)

dx9diff12=pu.get_surv_diff(nside=128)
dx9diff34=pu.get_surv_diff(surv1=3,surv2=4,nside=128)

dx91234= hp.ma((dx9diff34+dx9diff12)/2.)
f=open('fsl_sim_diffs_all_128.pkl','rb')
allsims=cPickle.load(f)
vallsims=np.array(allsims).T
f.close()
fsld27m0=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=0)
fsld27s0=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=0)
fsld28m0=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=0)
fsld28s0=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=0)

d27m0=hp.ma(hp.ud_grade(maskmap(fsld27m0),128,pess=True))
d27s0=hp.ma(hp.ud_grade(maskmap(fsld27s0),128,pess=True))
d28m0=hp.ma(hp.ud_grade(maskmap(fsld28m0),128,pess=True))
d28s0=hp.ma(hp.ud_grade(maskmap(fsld28s0),128,pess=True))

blk=hp.ma(range(196608))
theta,phi=hp.pix2ang(128,blk)
b=180.*(theta-np.pi/2.)/np.pi
blk=hp.ma(np.ones(196608,dtype=float))
blk.mask=np.abs(b)<10
totalmask128=dx91234.mask|blk.mask|d28m0.mask|d28s0.mask|d27s0.mask|d27m0.mask
cutvallsims=vallsims[totalmask128==False,:]
l27m1234=(l27m12+l27m34)/2.0
l27s1234=(l27s12+l27s34)/2.0
l28m1234=(l28m12+l28m34)/2.0
l28s1234=(l28s12+l28s34)/2.0
l27ms1234=(l27m1234+l27s1234)/2.0
l28ms1234=(l28m1234+l28s1234)/2.0

hp.mollview(l27ms1234,title='DX9 30 GHz 27MS 1234',min=-20e-6,max=10e-6)
hp.mollview(l28ms1234,title='DX9 30 GHz 28MS 1234',min=-20e-6,max=10e-6)

l27ms1234cut=np.array(l27ms1234[totalmask128==False])
l28ms1234cut=np.array(l28ms1234[totalmask128==False])

params27ms=np.linalg.lstsq(cutvallsims[:,0:7],l27ms1234cut)
params28ms=np.linalg.lstsq(cutvallsims[:,14:21],l28ms1234cut)
params27ms[0]
params28ms[0]

model27=np.sum(vallsims[:,0:7]*params27ms[0],axis=1)
model28=np.sum(vallsims[:,14:21]*params28ms[0],axis=1)

hp.mollview(model27,min=-20e-6,max=5e-6,title='Model 27 MS')
hp.mollview(model28,min=-20e-6,max=5e-6,title='Model 28 MS')

hp.mollview(l27ms1234-model27,title='DX9 30 GHz 27MS 1234-model27',min=-20e-6,max=5e-6)
hp.mollview(l28ms1234-model28,title='DX9 30 GHz 28MS 1234-model28',min=-20e-6,max=5e-6)

import cPickle
import healpy as hp
import planck_util as pu
dx9diff12=pu.get_surv_diff(nside=128)
dx9diff34=pu.get_surv_diff(surv1=3,surv2=4,nside=128)
dx91234= hp.ma((dx9diff34+dx9diff)/2.)
hp.mollview(dx91234,title='DX9 30 GHz mean (ss1-ss2)/2 , (ss3-ss4)/2',min=-20e-6,max=5e-6)
f=open('fsl_sim_diffs_all_128.pkl','rb')
allsims=cPickle.load(f)
vallsims=np.array(allsims).T

fsld27m0=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=0)
fsld27s0=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=0)
fsld28s0=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=0)
fsld28m0=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=0)

def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()
d28m0=hp.ma(hp.ud_grade(maskmap(fsld28m0),128,pess=True))
d28s0=hp.ma(hp.ud_grade(maskmap(fsld28s0),128,pess=True))
d27s0=hp.ma(hp.ud_grade(maskmap(fsld27s0),128,pess=True))
d27m0=hp.ma(hp.ud_grade(maskmap(fsld27m0),128,pess=True))

blk=hp.ma(range(196608))
theta,phi=hp.pix2ang(128,blk)
b=180.*(theta-np.pi/2.)/np.pi
blk=hp.ma(np.ones(196608,dtype=float))
blk.mask=np.abs(b)<10
totalmask128=dx91234.mask|blk.mask|d28m0.mask|d28s0.mask|d27s0.mask|d27m0.mask
cutvallsims=vallsims[totalmask128==False,:]
dx91234cut=np.array(dx91234[totalmask128==False])
ptest=np.linalg.lstsq(cutvallsims,dx91234cut)
type(ptest)
ptest[0]
model=np.sum(vallsims*ptest[0],axis=1)
hp.mollview(model,min=-20e-6,max=5e-6)
hp.mollview(model,min=-20e-6,max=5e-6,title='Model')
hp.mollview(dx91234-model,title='DX9 30 GHz mean (ss1-ss2)/2 , (ss3-ss4)/2  - Model',min=-20e-6,max=5e-6)

model.mask=totalmask128
hp.mollview(dx91234-model,title='DX9 30 GHz mean (ss1-ss2)/2 , (ss3-ss4)/2  - Model',min=-20e-6,max=5e-6)

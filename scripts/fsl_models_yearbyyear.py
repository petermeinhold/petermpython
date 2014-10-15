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

dx912= hp.ma(dx9diff12)
dx934= hp.ma(dx9diff34)
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
totalmask128a=dx912.mask|blk.mask|d28m0.mask|d28s0.mask|d27s0.mask|d27m0.mask
totalmask128b=dx934.mask|blk.mask|d28m0.mask|d28s0.mask|d27s0.mask|d27m0.mask

cutvallsimsa=vallsims[totalmask128a==False,:]
cutvallsimsb=vallsims[totalmask128b==False,:]

l27ms12=(l27m12+l27s12)/2.0
l28ms12=(l28m12+l28s12)/2.0
l27ms34=(l27m34+l27s34)/2.0
l28ms34=(l28m34+l28s34)/2.0

hp.mollview(l27ms12,title='DX9 30 GHz 27MS 12',min=-10e-6,max=10e-6)
hp.mollview(l28ms12,title='DX9 30 GHz 28MS 12',min=-10e-6,max=10e-6)

hp.mollview(l27ms34,title='DX9 30 GHz 27MS 34',min=-10e-6,max=10e-6)
hp.mollview(l28ms34,title='DX9 30 GHz 28MS 34',min=-10e-6,max=10e-6)


l27ms12cut=np.array(l27ms12[totalmask128a==False])
l28ms12cut=np.array(l28ms12[totalmask128a==False])

l27ms34cut=np.array(l27ms34[totalmask128b==False])
l28ms34cut=np.array(l28ms34[totalmask128b==False])

params27msa=np.linalg.lstsq(cutvallsimsa[:,0:7],l27ms12cut)
params28msa=np.linalg.lstsq(cutvallsimsa[:,14:21],l28ms12cut)
params27msb=np.linalg.lstsq(cutvallsimsb[:,0:7],l27ms34cut)
params28msb=np.linalg.lstsq(cutvallsimsb[:,14:21],l28ms34cut)

model27a=np.sum(vallsims[:,0:7]*params27msa[0],axis=1)
model28a=np.sum(vallsims[:,14:21]*params28msa[0],axis=1)

model27b=np.sum(vallsims[:,0:7]*params27msb[0],axis=1)
model28b=np.sum(vallsims[:,14:21]*params28msb[0],axis=1)

hp.mollview(model27a,min=-10e-6,max=10e-6,title='Model 27a MS')
hp.mollview(model28a,min=-10e-6,max=10e-6,title='Model 28a MS')
hp.mollview(model27b,min=-10e-6,max=10e-6,title='Model 27b MS')
hp.mollview(model28b,min=-10e-6,max=10e-6,title='Model 28b MS')

hp.mollview(l27ms12-model27a,title='DX9 30 GHz 27MS 12-model27',min=-10e-6,max=10e-6)
hp.mollview(l28ms12-model28a,title='DX9 30 GHz 28MS 12-model28',min=-10e-6,max=10e-6)
hp.mollview(l27ms34-model27b,title='DX9 30 GHz 27MS 34-model27',min=-10e-6,max=10e-6)
hp.mollview(l28ms34-model28b,title='DX9 30 GHz 28MS 34-model28',min=-10e-6,max=10e-6)

f=open('../frequency_bandpasses/bandpassweights_2.csv','r')
b=f.readlines()
f.close()
bb=np.zeros([7,5],dtype=float)
for i in range(7):
    bb[i,:]=np.array(np.str.split(b[i+1],','))
bbb=bb[:,1:].T
bbb=bbb.reshape(28)
bbb27=np.zeros(28)
bbb28=np.zeros(28)
bbb27[0:14]=bbb[0:14]
bbb28[14:]=bbb[14:]

bp27model=np.sum(vallsims*bbb27,axis=1)
bp28model=np.sum(vallsims*bbb28,axis=1)

hp.mollview(l27ms12-bp27model,title='DX9 30 GHz 27MS 12-bandpass model27',min=-10e-6,max=10e-6)
hp.mollview(l28ms12-bp28model,title='DX9 30 GHz 28MS 12-bandpass model28',min=-10e-6,max=10e-6)
hp.mollview(l27ms34-bp27model,title='DX9 30 GHz 27MS 34-bandpass model27',min=-10e-6,max=10e-6)
hp.mollview(l28ms34-bp28model,title='DX9 30 GHz 28MS 34-bandpass model28',min=-10e-6,max=10e-6)


hp.mollview(l27ms12,title='DX9 30 GHz 27MS 12',min=-10e-6,max=10e-6)
hp.mollview(l28ms12,title='DX9 30 GHz 28MS 12',min=-10e-6,max=10e-6)
hp.mollview(l27ms34,title='DX9 30 GHz 27MS 34',min=-10e-6,max=10e-6)
hp.mollview(l28ms34,title='DX9 30 GHz 28MS 34',min=-10e-6,max=10e-6)



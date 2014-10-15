# IPython log file

from petermpython import planck_util as pu

from glob import glob
import healpy as hp

nside=512

base30fsl='/global/scratch/sd/planck/user/peterm/sidelobes/sims_30ghz_bandpass/average_beam/'

fl30ave=glob(base30fsl+'*.fits')
fl30dip=glob(base30fsl+'dipole/*.fits')
fsl30gals1t=hp.ma(hp.read_map(fl30ave[1],0))
fsl30gals2t=hp.ma(hp.read_map(fl30ave[0],0))
fsl30gals1q=hp.ma(hp.read_map(fl30ave[1],1))
fsl30gals2q=hp.ma(hp.read_map(fl30ave[0],1))
fsl30gals1u=hp.ma(hp.read_map(fl30ave[1],2))
fsl30gals2u=hp.ma(hp.read_map(fl30ave[0],2))

fsl30dips1t=hp.ma(hp.read_map(fl30dip[0],0))
fsl30dips2t=hp.ma(hp.read_map(fl30dip[1],0))
fsl30dips1q=hp.ma(hp.read_map(fl30dip[0],1))
fsl30dips2q=hp.ma(hp.read_map(fl30dip[1],1))
fsl30dips1u=hp.ma(hp.read_map(fl30dip[0],2))
fsl30dips2u=hp.ma(hp.read_map(fl30dip[1],2))



#mask=dx9_30_ss12difft.mask | dx9_30_ss12diffq.mask | dx9_30_ss12diffu.mask |(abs(fsl30s2q)>1)|(abs(fsl30s1q)>1) 

#fsl30gals1t.mask=mask
#fsl30gals1q.mask=mask
#fsl30gals1u.mask=mask
#fsl30gals1t.mask=mask
#fsl30gals1q.mask=mask
#fsl30gals1u.mask=mask
#fsl30dips1t.mask=mask
#fsl30dips1q.mask=mask
#fsl30dips1u.mask=mask
#fsl30dips1t.mask=mask
#fsl30dips1q.mask=mask
#fsl30dips1u.mask=mask

fsl30s1t=fsl30gals1t+fsl30dips1t
fsl30s1q=fsl30gals1q+fsl30dips1q
fsl30s1u=fsl30gals1u+fsl30dips1u
fsl30s2t=fsl30gals2t+fsl30dips2t
fsl30s2q=fsl30gals2q+fsl30dips2q
fsl30s2u=fsl30gals2u+fsl30dips2u
fsl30t=(fsl30s1t+fsl30s2t)/2.
fsl30q=(fsl30s1q+fsl30s2q)/2.
fsl30u=(fsl30s1u+fsl30s2u)/2.
tmask=fsl30s1q.mask|fsl30s2q.mask 
fsl30t.mask=tmask
fsl30q.mask=tmask
fsl30u.mask=tmask

frac=float(len(fsl30q[fsl30q.mask==False]))/len(fsl30q)

#hp.mollview(fsl30t,min=-10e-6,max=10e-6,title='30 GHz band averaged FSL model, Survey 1+2 T')
#hp.mollview(fsl30q,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1+2 Q')
#hp.mollview(fsl30u,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1+2 U')

fsl30tqu=[fsl30t,fsl30q,fsl30u]
clfsltt,clfslee,clfslbb,clfslte,clfsltb,clfsleb=hp.anafast(fsl30tqu)
figure()
l=arange(len(clfsltt))
plot l,l*(l+1)*clfsltt/frac,label='FSL TT'
plot l,l*(l+1)*clfslee/frac,label='FSL EE'
plot l,l*(l+1)*clfslbb/frac,label='FSL BB'

title('FSL 30 GHz Survey 1+2 FSL sims'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2') 

fsl30galt=(fsl30gals1t+fsl30gals2t)/2.0
fsl30galq=(fsl30gals1q+fsl30gals2q)/2.0
fsl30galu=(fsl30gals1u+fsl30gals2u)/2.0
fsl30galt.mask=tmask
fsl30galq.mask=tmask
fsl30galu.mask=tmask
fsl30gal=[fsl30galt,fsl30galq,fsl30galu]

clfslgaltt,clfslgalee,clfslgalbb,clfslgalte,clfslgaltb,clfslgaleb=hp.anafast(fsl30gal)
figure()
plot l,l*(l+1)*clfslgaltt/frac,label='FSL gal TT'
plot l,l*(l+1)*clfslgalee/frac,label='FSL gal EE'
plot l,l*(l+1)*clfslgalbb/frac,label='FSL gal BB'

title('FSL 30 GHz Survey 1+2 FSL Gal only sims'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')
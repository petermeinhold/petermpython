#mess about to estimat N spectrum of WMAP 9 year
fwmapspec='/planck/WMAP\\wmap_tt_spectrum_9yr_v5.txt'
from glob import glob
import cPickle
import healpy as hp
import pyfits
import planck_util as pu
#aggregate sensitivities in mKsqrt(sec)
wvsens=1.13/sqrt(4)
wwsens=1.48/sqrt(8)
pvsens=.152
pspec=genfromtxt('wmap/base_planck_lowl_lowLike_highL_post_lensing.bestfit_cl')
fnom70='wmap/LFI_SkyMap_070_1024_R1.10_nominal.fits'
m70=pyfits.open(fnom70)
hits70=m70[1].data['Hits']
hits70=hp.ud_grade(hits70,nside_out=512,order_in='Nest',order_out='Ring')
planckweightmap=hits70/np.mean(hits70)

fwmv='wmap/wmap_band_imap_r9_9yr_V_v5.fits'
wmapw=pyfits.open(fwmv)

fwmap='wmap/wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat'
wmapcl=genfromtxt(fwmap)

wmapweightmap=wmapw[1].data['N_OBS']/np.mean( wmapw[1].data['N_OBS'])
#weighted time per pixel for  9year wmap
tpp_wmap_weighted=wmapweightmap*9*365.25*24*3600./hp.nside2npix(512)

#unweighted planck time per pixel
tpp_planck_nominal=planckweightmap*(15./12.)*365.25*24*3600./hp.nside2npix(512)
tpp_planck_mission=planckweightmap*4*365.25*24*3600./hp.nside2npix(512)
clsigwv_weighted=[]
clsigpvnom=[]
clsigpvmiss=[]
for i in range(2):
    mcw=1000.*randn(hp.nside2npix(512))*wwsens/sqrt(tpp_wmap_weighted)
    mcv=1000.*randn(hp.nside2npix(512))*wvsens/sqrt(tpp_wmap_weighted)
    mcwv=(mcw+mcv)/2.0
    mcpvnom=1000*randn(hp.nside2npix(512))*pvsens/sqrt(tpp_planck_nominal)
    mcpvmiss=1000*randn(hp.nside2npix(512))*pvsens/sqrt(tpp_planck_mission)
    clsigwv_weighted.append(hp.anafast(mcwv))
    clsigpvnom.append(hp.anafast(mcpvnom))
    clsigpvmiss.append(hp.anafast(mcpvmiss))
clsigwv_mean=np.mean(clsigwv_weighted,axis=0)
clsigpvnom_mean=np.mean(clsigpvnom,axis=0)
clsigpvmiss_mean=np.mean(clsigpvmiss,axis=0)

f70=open('ddx9_70ghz_halfdiff_cl_allsky.pkl','rb')
cl70hrhd=cPickle.load(f70)

l=arange(1536)+2

plot(l,l*(l+1)*clsigwv_mean/(2*np.pi),label='WMAP w+v, 9 yr, hit weighted wn average')
plot(l,l*(l+1)*clsigpvnom_mean/(2*np.pi),label='Planck 70 GHz, nominal, hit weighted wn average')
plot(l,l*(l+1)*clsigpvmiss_mean/(2*np.pi),label='Planck 70 GHz, 8 surveys, hit weighted wn average')
plot( wmapcl[:,0],wmapcl[:,2],label='WMAP noise (w+v) 9yr from likelihood, (EK)')
xlabel('Multipole',fontsize=18),ylabel('$l(l+1)C_l/2\pi , \muK^2$',fontsize=18),title('Simple white noise map comparisons',fontsize=20)
leg=legend(loc=2)
grid()

#old calculation, now have WMAP noise bias directly from Eiichiro:
#Attached please find the power spectrum file used in the WMAP 9-year likelihood code.
# (wmap/wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat)
#1st: multipole
#2nd: temperature power spectrum
#3rd: noise power spectrum
#4th: fsky

#sn=wspec[:,3]
#l1=arange(2,501)
#fsky1 = 0.826 - 0.091*(l1/500.)**2
#l2=arange(501,1201)
#fsky2 = 0.777 - 0.127*(500./l2)
#fsky=np.concatenate([fsky1,fsky2])
#l=np.concatenate([l1,l2])
#cth=wspec[:,1]*(2*np.pi/(l*(l+1)))
#cth=pspec[:1199,1]*(2*np.pi/(l*(l+1)))
#pi=np.pi
#nl=.25*(np.sqrt((4*cth)**2+4*2*(2*l+1)*((2*fsky*sn*pi/(l*(l+1)))**2))-4*cth)
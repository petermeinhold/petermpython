# IPython log file

from petermpython import planck_util as pu

from glob import glob
import healpy as hp
dx9_30_ss12difft=pu.get_surv_diff(nside=512)
dx9_30_ss12diffq=pu.get_surv_diff(nside=512,pol1='Q',pol2='Q')
dx9_30_ss12diffu=pu.get_surv_diff(nside=512,pol1='U',pol2='U')

nside=512
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)

base30fsl='/global/scratch/sd/planck/user/peterm/sidelobes/sims_30ghz_bandpass/average_beam/'

fl30ave=glob(base30fsl+'*.fits')
fl30dip=glob(base30fsl+'dipole/*.fits')
fsl30s1t=hp.ma(hp.read_map(fl30ave[1],0))
fsl30s2t=hp.ma(hp.read_map(fl30ave[0],0))
fsl30s1q=hp.ma(hp.read_map(fl30ave[1],1))
fsl30s1u=hp.ma(hp.read_map(fl30ave[1],2))
fsl30s2q=hp.ma(hp.read_map(fl30ave[0],1))
fsl30s2u=hp.ma(hp.read_map(fl30ave[0],2))

fsl30dips1t=hp.ma(hp.read_map(fl30dip[0],0))
fsl30dips2t=hp.ma(hp.read_map(fl30dip[1],0))
fsl30dips1q=hp.ma(hp.read_map(fl30dip[0],1))
fsl30dips1u=hp.ma(hp.read_map(fl30dip[0],2))
fsl30dips2u=hp.ma(hp.read_map(fl30dip[1],2))
fsl30dips2q=hp.ma(hp.read_map(fl30dip[1],1))

mask=dx9_30_ss12difft.mask | dx9_30_ss12diffq.mask | dx9_30_ss12diffu.mask |(abs(fsl30s2q)>1)|(abs(fsl30s1q)>1) 

fsl30gals1t.mask=mask
fsl30gals1q.mask=mask
fsl30gals1u.mask=mask
fsl30gals1t.mask=mask
fsl30gals1q.mask=mask
fsl30gals1u.mask=mask
fsl30dips1t.mask=mask
fsl30dips1q.mask=mask
fsl30dips1u.mask=mask
fsl30dips1t.mask=mask
fsl30dips1q.mask=mask
fsl30dips1u.mask=mask

fsl30s1t=fsl30gals1t+fsl30dips1t
fsl30s1q=fsl30gals1q+fsl30dips1q
fsl30s1u=fsl30gals1u+fsl30dips1u
fsl30s2t=fsl30gals2t+fsl30dips2t
fsl30s2q=fsl30gals2q+fsl30dips2q
fsl30s2u=fsl30gals2u+fsl30dips2u
fsl30t=fsl30s1t+fsl30s2t
fsl30q=fsl30s1q+fsl30s2q
fsl30u=fsl30s1u+fsl30s2u

fsl30difft=fsl30s1t-fsl30s2t
fsl30diffq=fsl30s1q-fsl30s2q
fsl30diffu=fsl30s1u-fsl30s2u

hp.mollview(fsl30s1t,min=-10e-6,max=10e-6,title='30 GHz band averaged FSL model, Survey 1  T')
show()
savefig('fsl_plots_912/fsl_ss1_t.png')
hp.mollview(fsl30s2t,min=-10e-6,max=10e-6,title='30 GHz band averaged FSL model, Survey 2  T')
show()
savefig('fsl_plots_912/fsl_ss2_t.png')
hp.mollview(fsl30s1q,min=-.1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1  Q')
show()
savefig('fsl_plots_912/fsl_ss1_q.png')
hp.mollview(fsl30s2q,min=-.1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 2  Q')
show()
savefig('fsl_plots_912/fsl_ss2_q.png')
hp.mollview(fsl30s1u,min=-.1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1  U')
show()
savefig('fsl_plots_912/fsl_ss1_u.png')
hp.mollview(fsl30s2u,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 2  U')
show()
savefig('fsl_plots_912/fsl_ss2_u.png')


hp.mollview(fsl30t,min=-10e-6,max=10e-6,title='30 GHz band averaged FSL model, Survey 1+2 T')
show()
savefig('fsl_plots_912/fsl_t.png')

hp.mollview(fsl30q,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1+2 Q')
show()
savefig('fsl_plots_912/fsl_q.png')

hp.mollview(fsl30u,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1+2 U')
show()
savefig('fsl_plots_912/fsl_u.png')

hp.mollview(fsl30difft,min=-10e-6,max=10e-6,title='30 GHz band averaged FSL model, Survey 1-2 T')
show()
savefig('fsl_plots_912/fsl_diff12_t.png')

hp.mollview(fsl30diffq,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1-2 Q')
show()
savefig('fsl_plots_912/fsl_diff12_q.png')

hp.mollview(fsl30diffu,min=-1e-6,max=1e-6,title='30 GHz band averaged FSL model, Survey 1-2 U')
show()
savefig('fsl_plots_912/fsl_diff12_u.png')

hp.mollview(dx9_30_ss12difft,min=-10e-6,max=10e-6,title='DX9 30 GHz Survey 1-2 T')
show()
savefig('dx9_30_ss12difft.png')
hp.mollview(dx9_30_ss12diffq,min=-10e-6,max=10e-6,title='DX9 30 GHz Survey 1-2 Q')
show()
savefig('dx9_30_ss12diffq.png')
hp.mollview(dx9_30_ss12diffu,min=-10e-6,max=10e-6,title='DX9 30 GHz Survey 1-2 U')
show()
savefig('dx9_30_ss12diffu.png')

dx9diff12_sub_fsldiff12t=dx9_30_ss12difft-fsl30difft
dx9diff12_sub_fsldiff12q=dx9_30_ss12diffu-fsl30diffq
dx9diff12_sub_fsldiff12u=dx9_30_ss12diffq-fsl30diffu

hp.mollview(dx9diff12_sub_fsldiff12t-mean(dx9diff12_sub_fsldiff12t[abs(dx9diff12_sub_fsldiff12t)<1]),min=-20e-6,max=20e-6,title='30 GHz Survey 1-2 T: DX9 - FSL sim')
show()
savefig('dx9_diff12_sub_fsl_diff12.png')

hp.mollview(dx9_30_ss12difft-mean(dx9_30_ss12difft[abs(dx9_30_ss12difft)<1]),min=-20e-6,max=20e-6,title='30 GHz Survey 1-2 T: DX9 ')
show()
savefig('dx9_diff12_sub_fsl_diff12.png')

dx9_30_ss12difft.mask=mask
dx9diff12_sub_fsldiff12t.mask=mask
dx9diff12_sub_fsldiff12q.mask=mask
dx9diff12_sub_fsldiff12u.mask=mask



specdx9diff12t=hp.anafast(dx9_30_ss12difft)
specdx9diff12t_subfslt=hp.anafast(dx9diff12_sub_fsldiff12t)

frac=len(dx9diff12_sub_fsldiff12t[dx9diff12_sub_fsldiff12t.mask==0])/float(hp.nside2npix(512))

l=arange(len(specdx9diff12t))
plot l,l*(l+1)*specdx9diff12t/frac,label='DX9 SS1-SS2'
plot l,l*(l+1)*specdx9diff12t_subfslt/frac,label='DX9-FSL sim SS1-SS2'
title('DX9 30 GHz Survey 1-2 TT spectrum vs FSL subtracted survey difference'),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2') 

dx9diff12tqu=[dx9_30_ss12difft,dx9_30_ss12diffq,dx9_30_ss12difft]
dx9diff12subfsltqu=[dx9diff12_sub_fsldiff12t,dx9diff12_sub_fsldiff12q,dx9diff12_sub_fsldiff12u]
cldx9ss12subfsltt,cldx9ss12subfslee,cldx9ss12subfslbb,cldx9ss12subfslte,cldx9ss12subfsltb,cldx9ss12subfsleb=hp.anafast(dx9diff12subfsltqu)
cldx9ss12tt,cldx9ss12ee,cldx9ss12bb,cldx9ss12te,cldx9ss12tb,cldx9ss12eb=hp.anafast(dx9diff12tqu)


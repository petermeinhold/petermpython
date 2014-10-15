# IPython log file

from petermpython import planck_util as pu
import healpy as hp
from glob import glob
freq=70
nside=1024
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
tmask=wmask.astype(np.bool)|psmask
frac=(1-float(tmask.sum())/len(tmask))

base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
fl30nominal=glob(base+'/*30*nominal*.fits')
basec='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/bandpass_correction/'
fl30nominalcorr=glob(basec+'/*correction*30*nominal*.fits')
map30nomtqu=hp.read_map(fl30nominal[0],[0,1,2])
corr30nomtqu=hp.read_map(fl30nominalcorr[0],[0,1,2])
corrmap30nom=[map30nomtqu[0]+corr30nomtqu[0],map30nomtqu[1]+corr30nomtqu[1],map30nomtqu[2]+corr30nomtqu[2]]

hp.mollview(corrmap30nom[0])
basefsl='gpeterm/sidelobes/ddx9/'
fl30nominalfsl=glob(basefsl+'*30*nom*4.fits')

fsl30nomtqu=hp.read_map(fl30nominalfsl[0],[0,1,2])
hp.mollview(fsl30nomtqu[0],title='FSL sim using DDX930 Ghz for input map, I'),min=-1e-6,max=1e-6
hp.mollview(fsl30nomtqu[1],title='FSL sim using DDX930 Ghz for input map, Q'),min=-1e-6,max=1e-6
hp.mollview(fsl30nomtqu[2],title='FSL sim using DDX930 Ghz for input map, u'),min=-1e-6,max=1e-6
fsl30nomtqumasked=hp.ma(fsl30nomtqu)
for map in fsl30nomtqumasked:
    map.mask=tmask

cl_ddx9_nom_fsl=hp.anafast(fsl30nomtqu)
cl_ddx9_nom_fsl_masked=hp.anafast(fsl30nomtqumasked)
l=arange(3072)
figure()

plot l,1e12*cl_ddx9_nom_fsl[0]*l*(l+1)/(2*np.pi),'b-',label='TT'
plot l,1e12*cl_ddx9_nom_fsl_masked[0]*l*(l+1)/(frac*2*np.pi),'b--',label='TT masked'

plot l,1e12*cl_ddx9_nom_fsl[1]*l*(l+1)/(2*np.pi),'g-',label='EE'
plot l,1e12*cl_ddx9_nom_fsl_masked[2]*l*(l+1)/(frac*2*np.pi),'g--',label='EE masked'

plot l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi),'r-',label='BB'
plot l,1e12*cl_ddx9_nom_fsl_masked[2]*l*(l+1)/(frac*2*np.pi),'r--',label='BB masked'

xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2') 
leg=legend()
pu.thicklegendlines(leg)
grid()
title('30 GHz bandweighted FSL model, Delta DX9 as input sky')


map30nomtqu_half1=hp.read_map(fl30nominal[1],[0,1,2])
map30nomtqu_half2=hp.read_map(fl30nominal[2],[0,1,2])
corrmap30nomtqu_half1=[map30nomtqu_half1[0]+corr30nomtqu[0],map30nomtqu_half1[1]+corr30nomtqu[1],map30nomtqu_half1[2]+corr30nomtqu[2]]
corrmap30nomtqu_half2=[map30nomtqu_half2[0]+corr30nomtqu[0],map30nomtqu_half2[1]+corr30nomtqu[1],map30nomtqu_half2[2]+corr30nomtqu[2]]

corrmap30nomtqu_half1=hp.ma(corrmap30nomtqu_half1)
corrmap30nomtqu_half2=hp.ma(corrmap30nomtqu_half2)

emask=np.where(corrmap30nomtqu_half1[0]==hp.UNSEEN)

tmask[emask[0]]=True

emask1=np.where(corrmap30nomtqu_half1[0]==hp.UNSEEN)
emask2=np.where(corrmap30nomtqu_half2[0]==hp.UNSEEN)
emask1q=np.where(corrmap30nomtqu_half1[1]==hp.UNSEEN)
emask1u=np.where(corrmap30nomtqu_half1[2]==hp.UNSEEN)
emask2u=np.where(corrmap30nomtqu_half2[2]==hp.UNSEEN)
emask2q=np.where(corrmap30nomtqu_half2[1]==hp.UNSEEN)
tmask[emask2[0]]=True
tmask[emask2q[0]]=True
tmask[emask2u[0]]=True
tmask[emask1u[0]]=True
tmask[emask1q[0]]=True
for map in corrmap30nomtqu_half1:
    map.mask=tmask
for map in corrmap30nomtqu_half2:
    map.mask=tmask
diff=[]
for map1,map2 in zip(corrmap30nomtqu_half1,corrmap30nomtqu_half2):
diff.append(map1-map2)
frac=(1-float(tmask.sum())/len(tmask))
cl_ddx9_30ghz_half1subhalf2=hp.anafast(diff)  
cl_ddx9_30ghz_half1xhalf2=hp.anafast(corrmap30nomtqu_half1,corrmap30nomtqu_half2)

figure()
plot l,1e12*cl_ddx9_nom_fsl[1]*l*(l+1)/(2*np.pi),'g-',label='FSL EE'
plot l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi),'r-',label='FSL BB'
plot l,1e12*cl_ddx9_30ghz_half1xhalf2[1]*l*(l+1)/(frac*2*np.pi),'c',label='DDX9 bpcorrected Half1x2 EE'
plot l,1e12*cl_ddx9_30ghz_half1xhalf2[2]*l*(l+1)/(frac*2*np.pi),'m',label='DDX9 bpcorrected Half1x2 BB'
plot l,1e12*cl_ddx9_30ghz_half1subhalf2[1]*l*(l+1)/(frac*2*np.pi),'b',label='DDX9 null Half1-half2 EE'
plot l,1e12*cl_ddx9_30ghz_half1subhalf2[2]*l*(l+1)/(frac*2*np.pi),'y',label='DDX9 null Half1-half2 BB'

xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2') 
leg=legend(loc=7)
pu.thicklegendlines(leg)
grid()
title('30 GHz bandweighted FSL model, Delta DX9 as input sky, compare with data')
xlim(0,50)
ylim(.0001,100)
yscale('log')
savefig('cl_30ghz_fsl_3_eebb_withdata_log.png')





fl30yr=glob(base+'/*30*yr*.fits')
fl30yr

fl70yr=glob(base+'/*70*yr*.fits')

fl70yr
map70y1tqu=hp.read_map(fl70yr[1],[0,1,2])
map70y2tqu=hp.read_map(fl70yr[0],[0,1,2])
map70y1tqu=hp.ma(map70y1tqu)
map70y2tqu=hp.ma(map70y2tqu)
for map in map70y1tqu:
    map.mask=tmask    
for map in map70y2tqu:
    map.mask=tmask







plot(l,cl_ddx9_nom_fsl[1]*l*(l+1)/(2*np.pi))
plot(l,cl_ddx9_nom_fsl[0]*l*(l+1)/(2*np.pi))
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi))
title='DDx9 30 GHz nominal FSL simulation')
xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2') 

xscale('log'),yscale('log')
figure(5)
xscale('log'),yscale('log')
xscale('linear'),yscale('linear')
figure()
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi))
autocalll 2
get_ipython().magic(u'autocall 2')
xscale('log'),yscale('log')
figure()
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'),xrange([0,12]))
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'),xrange=[0,12])
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'),xrange(0,12))
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'))
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'))
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'))
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'))
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi))
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi),xlabel('l'))
hold(False)
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi))
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'))
from matplotlib reload xlabel
from matplotlib load xlabel
from matplotlib import xlabel
from plt.matplotlib import xlabel
from matplotlib import xlabel
get_ipython().magic(u'whos ')
plt.xlabel()
get_ipython().magic(u'pinfo plt.xlabel')
get_ipython().magic(u'pinfo plt.ylabel')
plt.ylabel()
get_ipython().magic(u'pinfo plt.ylabel')
plot(l,cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi),xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2'))
plot(l,1e12*cl_ddx9_nom_fsl[2]*l*(l+1)/(2.*np.pi))
xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2')
xrange(0,12)
xrange([0,12])
xrange(0,12)
show()
get_ipython().magic(u'pinfo xrange')
q
get_ipython().system(u'clear ')
xlim(0,12)
ylim(0,.5)
plot(l,1e12*cl_ddx9_nom_fsl[1]*l*(l+1)/(2.*np.pi))
xlim(0,12)
ylim(0,.5)
plot(l,1e12*cl_ddx9_nom_fsl[0]*l*(l+1)/(2.*np.pi))
xlim(0,12)
ylim(0,.5)
from petermpython import planck_util as pu
import healpy as hp
from glob import glob
base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
fla=glob(base+'LFI*19_22*survey*.fits')
fla.sort()
flafull=glob(base+'LFI*19_22_full*.fits')
dbs1=pu.get_surv_diff_fname(fla[0],flafull[0],fwhm=1.,freq=70)
dbs2=pu.get_surv_diff_fname(fla[1],flafull[0],fwhm=1.,freq=70)
dbs3=pu.get_surv_diff_fname(fla[2],flafull[0],fwhm=1.,freq=70)
dbs4=pu.get_surv_diff_fname(fla[3],flafull[0],fwhm=1.,freq=70)
dbs5=pu.get_surv_diff_fname(fla[4],flafull[0],fwhm=1.,freq=70)
fl70year=glob(base+'/*70*year*.fits')
fl70year
fl70year=glob(base+'/*70*y*.fits')
fl70year
fl70yr
fl70yr=glob(base+'/*70*yr*.fits')
fl70yr
map70y1tqu=hp.read_map(fl70yr[1],[0,1,2])
map70y2tqu=hp.read_map(fl70yr[0],[0,1,2])
hp.mollview(map70y1tqu[1])
cl70y1x2=hp.anafast(map70y1tqu,map70y2tqu)
cl70y1x1=hp.anafast(map70y1tqu)
cl70y2x2=hp.anafast(map70y2tqu)
len(cl_ddx9_nom_fsl[0])
len(cl70y2x2[0])
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY1 EE')
figure(6)
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY1 EE')
hold(True)
xlim(0,12)
ylim(0,.5)
ylim(0,5)
ylim(0,35)
ylim(0,15)
plot(l,1e12*cl70y2x2[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y2xY2 EE')
plot(l,1e12*cl70y1x2[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY2 EE')
get_ipython().magic(u'pinfo hp.anafast')
get_ipython().magic(u'pinfo pu.make_bmask')
wmask
get_ipython().set_next_input(u'bmask=pu.make_mask');get_ipython().magic(u'pinfo pu.make_mask')
get_ipython().magic(u'pinfo pu.make_mask')
get_ipython().magic(u'pinfo pu.make_bmask')
bmask=pu.make_bmask(bmax=20)
type(map70y2tqu[0])
for map in map70y2tqu:
    map=hp.ma(map)
    map.mask=bmask
    
for map in map70y1tqu:
    map=hp.ma(map)
    map.mask=bmask
    
cl70y1x2=hp.anafast(map70y1tqu,map70y2tqu)
cl70y1x1=hp.anafast(map70y1tqu)
cl70y2x2=hp.anafast(map70y2tqu)
figure()
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY1 EE')
ylim(0,15)
xlim(0,12)
plot(l,1e12*cl70y1x2[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY2 EE')
plot(l,1e12*cl70y2x2[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y2xY2 EE')
leg=legend()
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
freq=70
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
nside=1024
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
tmask=bmask|psmask
bmask.dtype
tmask=bmask astype np.bool|psmask
tmask=bmask astype np.boolean|psmask
tmask=bmask astype np.bool()|psmask
tmask=bmask astype( np.bool)|psmask
tmask=bmask.astype(np.bool)|psmask
hp.mollview(tmask)
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
type(wmask)
tmask=wmask.astype(np.bool)|psmask
hp.mollview(tmask)
for map in map70y1tqu:
    map=hp.ma(map)
    map.mask=tmask
    
for map in map70y2tqu:
    map=hp.ma(map)
    map.mask=tmask
    
cl70y1x2=hp.anafast(map70y1tqu,map70y2tqu)
type(map70y1tqu[0])
for map in map70y2tqu:
    map[tmask]=hp.UNSEEN
    
for map in map70y1tqu:
    map[tmask]=hp.UNSEEN
    
hp.mollview(map70y1tqu[0])
cl70y1x2=hp.anafast(map70y1tqu,map70y2tqu)
cl70y1x1=hp.anafast(map70y1tqu)
cl70y2x2=hp.anafast(map70y2tqu)
figure()
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY1 EE')
ylim(0,15)
xlim(0,12)
ylim(0,2)
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y2xY2 EE')
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY2 EE')
figure()
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY1 EE')
plot(l,1e12*cl70y2x2[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y2xY2 EE')
plot(l,1e12*cl70y1x2[1]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY2 EE')
ylim(0,2)
xlim(0,2)
xlim(0,12)
leg=legend()
grid()
ylim(0,.5)
ylim(0,1)
title('DDX9 70 GHz, Y1, Y2 EE, wmap and ptsrc masked, no correction for skyfrac')
xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2')
figure()
plot(l,1e12*cl70y1x1[2]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY1 BB')
plot(l,1e12*cl70y2x2[2]*l*(l+1)/(2.*np.pi),label='70 GHz Y2xY2 BB')
plot(l,1e12*cl70y1x2[2]*l*(l+1)/(2.*np.pi),label='70 GHz Y1xY2 BB')
ylim(0,1)
xlim(0,12)
xlim(0,20)
xlim(0,50)
y
ylim(-1,1)
grid()
ylim(-.5,.5)
xlim(0,12)
len(tmask==True)/len(tmask)
float(len(tmask==True))/len(tmask)
float(len(tmask.sum()))/len(tmask)
float(tmask.sum())/len(tmask)
1/(1-float(tmask.sum())/len(tmask))
frac=(1-float(tmask.sum())/len(tmask))
figure()
plot(l,1e12*cl70y1x2[2]*l*(l+1)/(frac*2.*np.pi),label='70 GHz Y1xY2 BB')
plot(l,1e12*cl70y1x1[2]*l*(l+1)/(frac*2.*np.pi),label='70 GHz Y1xY1 BB')
plot(l,1e12*cl70y2x2[2]*l*(l+1)/(frac*2.*np.pi),label='70 GHz Y2xY2 BB')
xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2')
title('DDX9 70 GHz, Y1, Y2 BB, wmap and ptsrc masked, pseudo-Cl corrected for skyfraction')
grid()
savefig('plots/oct8/ddx9_70ghz_bb_wide.png')
xlim(0,50)
ylim(-5,5)
savefig('plots/oct8/ddx9_70ghz_bb_med.png')
xlim(0,12)
ylim(-.5,.5)
ylim(-1,1)
ylim(-.5,1)
savefig('plots/oct8/ddx9_70ghz_bb_small.png')
figure()
figure(7)
leg=legend()
pu.thicklegendlines()
pu.thicklegendlines(leg)
savefig('plots/oct8/ddx9_70ghz_bb_small.png')
xlim(0,50)
ylim(-5,5)
savefig('plots/oct8/ddx9_70ghz_bb_med.png')
xlim(0,1500)
ylim(-50,50)
ylim(-500,500)
ylim(-5000,5000)
ylim(-500,5000)
savefig('plots/oct8/ddx9_70ghz_bb_wide.png')
figure(8)
plot(l,1e12*cl70y1x2[1]*l*(l+1)/(frac*2.*np.pi),label='70 GHz Y1xY2 EE')
plot(l,1e12*cl70y1x1[1]*l*(l+1)/(frac*2.*np.pi),label='70 GHz Y1xY1 EE')
plot(l,1e12*cl70y2x2[1]*l*(l+1)/(frac*2.*np.pi),label='70 GHz Y2xY2 EE')
leg=legend()
pu.thicklegendlines(leg)
title('DDX9 70 GHz, Y1, Y2 EE, wmap and ptsrc masked, pseudo-Cl corrected for skyfraction')
xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2')
grid()
xlim(0,1500)
ylim(-500,5000)
savefig('plots/oct8/ddx9_70ghz_ee_wide.png')
xlim(0,50)
ylim(-5,5)
savefig('plots/oct8/ddx9_70ghz_ee_med.png')
xlim(0,12)
ylim(-.5,1)
savefig('plots/oct8/ddx9_70ghz_ee_small.png')

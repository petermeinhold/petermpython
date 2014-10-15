import matplotlib.pyplot as plt
import healpy as hp
from glob import glob
import planck_util as pu
import cPickle
basedelta='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
#crap to create rectangular image to save to  pkl to print on windows!
nside=1024
xsize = 2000
ysize = xsize/2.
unit = "mK"
vmin = -1; vmax = 1
theta = np.linspace(np.pi, 0, ysize)
phi   = np.linspace(-np.pi, np.pi, xsize)
longitude = np.radians(np.linspace(-180, 180, xsize))
latitude = np.radians(np.linspace(-90, 90, ysize))
PHI, THETA = np.meshgrid(phi, theta)

#explicit filenames to avoid confusion
f30rh1= '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_1_nominal_1s.fits'
f30rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_ringhalf_2_nominal_1s.fits'
f30ss1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_survey_1_1s.fits'
f30ss2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_survey_2_1s.fits'
f30ss3='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_survey_3_1s.fits'

f44rh1= '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_1_nominal_1s.fits'
f44rh2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_ringhalf_2_nominal_1s.fits'
f44ss1='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_survey_1_1s.fits'
f44ss2='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_survey_2_1s.fits'
f44ss3='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_survey_3_1s.fits'

f70rh1 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_1_nominal_1s.fits'
f70rh2 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_2_nominal_1s.fits'
f70ss1 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_survey_1_1s.fits'
f70ss2 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_survey_2_1s.fits'
f70ss3 = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_survey_3_1s.fits'

#nside=1024
#psmask30 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_30GHz_*.*')[0]), nside)))
#psmask44 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_44GHz_*.*')[0]), nside)))
#psmask70 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/sources_lfi/mask_ps_70GHz_*.*')[0]), nside)))

#mask80 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/global/homes/p/peterm/masks/galactic_new/*.80_sky*.fits')[0]), nside)))
#tmask30=psmask30.astype(np.bool) #|mask80.astype(np.bool)
#tmask44=psmask44.astype(np.bool) #|mask80.astype(np.bool)
#tmask70=psmask70.astype(np.bool) #|mask80.astype(np.bool)

fwhm=2.0*np.pi/180.

m30rh1=hp.ma(hp.read_map(f30rh1,[0]))
m30rh2=hp.ma(hp.read_map(f30rh2,[0]))
m30ss1=hp.ma(hp.read_map(f30ss1,[0]))
m30ss2=hp.ma(hp.read_map(f30ss2,[0]))
m30ss3=hp.ma(hp.read_map(f30ss3,[0]))

d30rh=(m30rh1-m30rh2)/2.
d30ss12=(m30ss1-m30ss2)/2.
d30ss13=(m30ss1-m30ss3)/2.
m30rh1=0
m30rh2=0
m30ss1=0
m30ss2=0
m30ss3=0
d30rh=hp.smoothing(d30rh,fwhm=fwhm)
d30ss12=hp.smoothing(d30ss12,fwhm=fwhm)
d30ss13=hp.smoothing(d30ss13,fwhm=fwhm)

d30rh=d30rh-np.mean(d30rh)
d30ss12=d30ss12-np.mean(d30ss12)
d30ss13=d30ss13-np.mean(d30ss13)

f=open('/global/homes/p/peterm/plots/clplots/d30ss12.pkl','wb')
grid_map = d30ss12[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

f=open('/global/homes/p/peterm/plots/clplots/d30ss13.pkl','wb')
grid_map = d30ss13[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

f=open('/global/homes/p/peterm/plots/clplots/d30rh.pkl','wb')
grid_map = d30rh[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()
d30rh=0
d30ss12=0
d30ss13=0

m44rh1=hp.ma(hp.read_map(f44rh1,[0]))
m44rh2=hp.ma(hp.read_map(f44rh2,[0]))
m44ss1=hp.ma(hp.read_map(f44ss1,[0]))
m44ss2=hp.ma(hp.read_map(f44ss2,[0]))
m44ss3=hp.ma(hp.read_map(f44ss3,[0]))

d44rh=(m44rh1-m44rh2)/2.
d44ss12=(m44ss1-m44ss2)/2.
d44ss13=(m44ss1-m44ss3)/2.
m44rh1=0
m44rh2=0
m44ss1=0
m44ss2=0
m44ss3=0
d44rh=hp.smoothing(d44rh,fwhm=fwhm)
d44ss12=hp.smoothing(d44ss12,fwhm=fwhm)
d44ss13=hp.smoothing(d44ss13,fwhm=fwhm)

d44rh=d44rh-np.mean(d44rh)
d44ss12=d44ss12-np.mean(d44ss12)
d44ss13=d44ss13-np.mean(d44ss13)


f=open('/global/homes/p/peterm/plots/clplots/d44ss12.pkl','wb')
grid_map = d44ss12[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

f=open('/global/homes/p/peterm/plots/clplots/d44ss13.pkl','wb')
grid_map = d44ss13[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

f=open('/global/homes/p/peterm/plots/clplots/d44rh.pkl','wb')
grid_map = d44rh[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

m70rh1=hp.ma(hp.read_map(f70rh1,[0]))
m70rh2=hp.ma(hp.read_map(f70rh2,[0]))
m70ss1=hp.ma(hp.read_map(f70ss1,[0]))
m70ss2=hp.ma(hp.read_map(f70ss2,[0]))
m70ss3=hp.ma(hp.read_map(f70ss3,[0]))

d70rh=(m70rh1-m70rh2)/2.
d70ss12=(m70ss1-m70ss2)/2.
d70ss13=(m70ss1-m70ss3)/2.
m70rh1=0
m70rh2=0
m70ss1=0
m70ss2=0
m70ss3=0
d70rh=hp.smoothing(d70rh,fwhm=fwhm)
d70ss12=hp.smoothing(d70ss12,fwhm=fwhm)
d70ss13=hp.smoothing(d70ss13,fwhm=fwhm)

d70rh=d70rh-np.mean(d70rh)
d70ss12=d70ss12-np.mean(d70ss12)
d70ss13=d70ss13-np.mean(d70ss13)

# map on a matrix

f=open('/global/homes/p/peterm/plots/clplots/d70ss12.pkl','wb')
grid_map = d70ss12[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

f=open('/global/homes/p/peterm/plots/clplots/d70ss13.pkl','wb')
grid_map = d70ss13[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()

f=open('/global/homes/p/peterm/plots/clplots/d70rh.pkl','wb')
grid_map = d70rh[hp.ang2pix(nside, THETA, PHI)]
cPickle.dump(grid_map,f)
f.close()



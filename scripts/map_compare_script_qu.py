sys.path.append('/global/homes/p/peterm/petermpython/paperplots/python/scripts')
sys.path.append('/global/homes/p/peterm/petermpython/paperplots/data')
import cPickle
from glob import glob
from setup_matplotlib import *
from matplotlib.colors import ListedColormap
import healpy as hp
import planck_util as pu

def sum_wmap_maps(filelist):
    m=hp.read_map(filelist[0],(0,1,2,3,4))
    msumi=m[0]*m[4]
    msumq=m[1]*m[4]
    msumu=m[2]*m[4]
    msums=m[3]*m[4]
    msumn=m[4]
    for f in filelist[1:]:
        m=hp.read_map(f,(0,1,2,3,4))
        msumi+=m[0]*m[4]
        msumq+=m[1]*m[4]
        msumu+=m[2]*m[4]
        msums+=m[3]*m[4]
        msumn+=m[4]
    msumi=msumi/msumn
    msumq=msumq/msumn
    msumu=msumu/msumn
    msums=msums/msumn
    return msumi,msumq,msumu,msums,msumn
    

f70_dx9= ['/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120911_full_1s.fits']
f70_dx10=['/global/project/projectdirs/planck/data/mission/DPC_maps/dx10_pre/RadiometerCalibration/70/LFI_70_1024_20130503_RAD_full_1s.fits']
flw= ['wmap/wmap_band_iqumap_r9_9yr_Q_v5.fits','wmap/wmap_band_iqumap_r9_9yr_V_v5.fits', 'wmap/wmap_band_iqumap_r9_9yr_W_v5.fits']
flh=['/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_100_2048_20130109_nominal_noZodi.fits',
'/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_143_2048_20130109_nominal_noZodi.fits',
'/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_217_2048_20130109_nominal_noZodi.fits']



cmap = ListedColormap(np.loadtxt("petermpython/paperplots/data/Planck_Parchment_RGB.txt")/255.)

fgalmask40=glob('/global/u1/p/peterm/masks/galactic_new/comb*.4*.fits')[0]
galmask40 =(hp.ud_grade(hp.read_map(fgalmask40),512) < 1).astype(np.bool)
fgalmask70=glob('/global/u1/p/peterm/masks/galactic_new/comb*.7*.fits')[0]
galmask70 =(hp.ud_grade(hp.read_map(fgalmask70),512) < 1).astype(np.bool)
freq=70
psmask70 = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), 512)))

mq=hp.ma(hp.read_map(flw[0],(0,1,2)))
mv=hp.ma(hp.read_map(flw[1],(0,1,2)))
mw=hp.ma(hp.read_map(flw[2],(0,1,2)))

m100=hp.ma(hp.read_map(flh[0],(0,1,2)))
m100t=hp.ud_grade(m100[0],nside_out=512)*1000000.
m100q=hp.ud_grade(m100[1],nside_out=512)*1000000.
m100u=hp.ud_grade(m100[2],nside_out=512)*1000000.



m70dx9=hp.ma(hp.read_map(f70_dx9[0],(0,1,2)))
m70dx10=hp.ma(hp.read_map(f70_dx10[0],(0,1,2)))
m70dx9t=hp.ud_grade(m70dx9[0],nside_out=512)*1000000.
m70dx9q=hp.ud_grade(m70dx9[1],nside_out=512)*1000000.
m70dx9u=hp.ud_grade(m70dx9[2],nside_out=512)*1000000.
m70dx10t=hp.ud_grade(m70dx10[0],nside_out=512)*1000000.
m70dx10q=hp.ud_grade(m70dx10[1],nside_out=512)*1000000.
m70dx10u=hp.ud_grade(m70dx10[2],nside_out=512)*1000000.

mqt=mq[0]*1000.
mqq=mq[1]*1000.
mqu=mq[2]*1000.

mvt=mv[0]*1000.
mvq=mv[1]*1000.
mvu=mv[2]*1000.

mwt=mw[0]*1000.
mwq=mw[1]*1000.
mwu=mw[2]*1000.


#tmask=galmask40|psmask70
tmask=galmask70|psmask70

mqt.mask=tmask
mqq.mask=tmask
mqu.mask=tmask

mvt.mask=tmask
mvq.mask=tmask
mvu.mask=tmask

mwt.mask=tmask
mwq.mask=tmask
mwu.mask=tmask

m70dx9t.mask=tmask
m70dx9q.mask=tmask
m70dx9u.mask=tmask

m70dx10t.mask=tmask
m70dx10q.mask=tmask
m70dx10u.mask=tmask

d70dx9vq=m70dx9q-mvq
d70dx9vu=m70dx9u-mvu

d70dx10vq=m70dx10q-mvq
d70dx10vu=m70dx10u-mvu

d70dx9vqsm=hp.smoothing(d70dx9vq,fwhm=radians(10.))
d70dx10vqsm=hp.smoothing(d70dx10vq,fwhm=radians(10.))

d70dx9wq=m70dx9q-mwq
d70dx9wu=m70dx9u-mwu

d70dx10wq=m70dx10q-mwq
d70dx10wu=m70dx10u-mwu

d70dx9100q=m70dx9q-m100q

dqvq=mqq-mvq
dqvqsm=hp.smoothing(dqvq,fwhm=radians(10.))

dvwq=mvq-mwq
dvwqsm=hp.smoothing(dvwq,fwhm=radians(10.))

d70dx9100qsm=hp.smoothing(d70dx9100q,fwhm=radians(10.))

d70dx9wqsm=hp.smoothing(d70dx9wq,fwhm=radians(10.))
d70dx10wqsm=hp.smoothing(d70dx10wq,fwhm=radians(10.))

hp.mollview d70dx9vqsm,min=-10,max=10,title='LFI70 ddx9 -WMAP Vband,  Q smooth 10 degrees'
hp.mollview d70dx9wqsm,min=-10,max=10,title='LFI70 ddx9 -WMAP Wband,  Q smooth 10 degrees'
hp.mollview d70dx10vqsm,min=-10,max=10,title='LFI70 ddx10 -WMAP Vband,  Q smooth 10 degrees'
hp.mollview d70dx10wqsm,min=-10,max=10,title='LFI70 ddx10 -WMAP Wband,  Q smooth 10 degrees'
hp.mollview dvwqsm,min=-10,max=10,title='WMAP Vband-WMAP Wband,  Q smooth 10 degrees'
hp.mollview dqvqsm,min=-10,max=10,title='WMAP Qband-WMAP Vband,  Q smooth 10 degrees'
hp.mollview d70dx9100qsm,min=-10,max=10,title='LFI70 - HFI 100,  Q smooth 10 degrees'



v_14=hp.ma(sum_wmap_maps(flwv[:8]))
v_58=hp.ma(sum_wmap_maps(flwv[8:16]))
v_14[0].mask=tmask
v_14[1].mask=tmask
v_14[2].mask=tmask
v_14[3].mask=tmask
v_58[0].mask=tmask
v_58[1].mask=tmask
v_58[2].mask=tmask
v_58[3].mask=tmask


tmask=galmask70|psmask70
vdiffi=v_14[0]-v_58[0]
vdiffq=v_14[1]-v_58[1]
vdiffu=v_14[2]-v_58[2]
vdiffs=v_14[3]-v_58[3]
vdiffi.mask=tmask
vdiffq.mask=tmask
vdiffu.mask=tmask
vdiffs.mask=tmask
vdiffp=np.sqrt(vdiffq**2+vdiffu**2)
vdiffism=hp.smoothing(vdiffi,fwhm=np.radians(10))
vdiffqsm=hp.smoothing(vdiffq,fwhm=np.radians(10))
vdiffusm=hp.smoothing(vdiffu,fwhm=np.radians(10))
vdiffpsm=hp.smoothing(vdiffp,fwhm=np.radians(10))
hp.mollview vdiffpsm*1000.,title='WMAP V band y1..4 - y5..8, quad sum Q and U',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_vband_p.png')
hp.mollview vdiffqsm*1000.,title='WMAP V band y1..4 - y5..8,  Q ',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_vband_q.png')
hp.mollview vdiffusm*1000.,title='WMAP V band y1..4 - y5..8,  U',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_vband_u.png')
hp.mollview vdiffism*1000.,title='WMAP V band y1..4 - y5..8,  I',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_vband_i.png')

w_14=hp.ma(sum_wmap_maps(flww[:16]))
w_58=hp.ma(sum_wmap_maps(flww[16:32]))



wdiffi=w_14[0]-w_58[0]
wdiffq=w_14[1]-w_58[1]
wdiffu=w_14[2]-w_58[2]
wdiffs=w_14[3]-w_58[3]
wdiffi=hp.ma(wdiffi)
wdiffq=hp.ma(wdiffq)
wdiffu=hp.ma(wdiffu)
wdiffs=hp.ma(wdiffs)
wdiffi.mask=tmask
wdiffq.mask=tmask
wdiffu.mask=tmask
wdiffs.mask=tmask
wdiffp=np.sqrt(wdiffq**2+wdiffu**2)
wdiffism=hp.smoothing(wdiffi,fwhm=np.radians(10))
wdiffqsm=hp.smoothing(wdiffq,fwhm=np.radians(10))
wdiffusm=hp.smoothing(wdiffu,fwhm=np.radians(10))
wdiffpsm=hp.smoothing(wdiffp,fwhm=np.radians(10))
hp.mollview wdiffpsm*1000.,title='WMAP WV band y1..4 - y5..8, quad sum Q and U',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_wband_p.png')
hp.mollview wdiffqsm*1000.,title='WMAP W band y1..4 - y5..8,  Q ',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_wband_q.png')
hp.mollview wdiffusm*1000.,title='WMAP W band y1..4 - y5..8,  U',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_wband_u.png')
hp.mollview wdiffism*1000.,title='WMAP W band y1..4 - y5..8,  I',unit='$\mu K_{cmb}$'
savefig('wmap_yrdiff_wband_i.png')





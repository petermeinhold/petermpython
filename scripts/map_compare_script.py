sys.path.append('/global/homes/p/peterm/petermpython/paperplots/python/scripts')
sys.path.append('/global/homes/p/peterm/petermpython/paperplots/data')
import cPickle
from glob import glob
from setup_matplotlib import *
from matplotlib.colors import ListedColormap
import healpy as hp
cmap = ListedColormap(np.loadtxt("petermpython/paperplots/data/Planck_Parchment_RGB.txt")/255.)

flw=['wmap/wmap_band_imap_r9_9yr_V_v5.fits', 'wmap/wmap_band_imap_r9_9yr_W_v5.fits']

fll=['/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_44_1024_20120914_nominal_1s.fits',
'/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_30_1024_20120914_nominal_1s.fits',
'/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_nominal_1s.fits']

flh=['/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_100_2048_20130109_nominal_noZodi.fits',
'/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_143_2048_20130109_nominal_noZodi.fits',
'/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_217_2048_20130109_nominal_noZodi.fits']


fgalmask40=glob('/global/u1/p/peterm/masks/galactic_new/comb*.4*.fits')[0]
galmask40 =(hp.ud_grade(hp.read_map(fgalmask40),512) < 1).astype(np.bool)
freq=30
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), 512)))

dvw=mv-mw
d3044=m30-m44
d4470=m44-m70
d70100=m70-m100
d100143=m100-m143
d143217=m143-m217
d44100=m44-m100
d70143=m70-m143

dv44=mv-m44
dv70=mv-m70
dv100=mv-m100
dv143=mv=m143

dw44=mw-m44
dw70=mw-m70
dw100=mw-m100
dw143=mw-m143

ngpvec=hp.ang2vec(0,np.pi)
sgpvec=hp.ang2vec(np.pi,np.pi)
ngppix=hp.query_disc(nside=512,vec=ngpvec,radius=np.radians(1))
sgppix=hp.query_disc(nside=512,vec=sgpvec,radius=np.radians(1))

mv=hp.ud_grade(hp.ma(hp.read_map(flw[0])),512,order_out='Ring')*1000.
mw=hp.ud_grade(hp.ma(hp.read_map(flw[1])),512,order_out='Ring')*1000.
m30=hp.ud_grade(hp.ma(hp.read_map(fll[0])),512)*1.e6
m44=hp.ud_grade(hp.ma(hp.read_map(fll[1])),512)*1.e6
m70=hp.ud_grade(hp.ma(hp.read_map(fll[2])),512)*1.e6
m100=hp.ud_grade(hp.ma(hp.read_map(flh[0])),512)*1.e6
m143=hp.ud_grade(hp.ma(hp.read_map(flh[1])),512)*1.e6
m217=hp.ud_grade(hp.ma(hp.read_map(flh[2])),512)*1.e6
mv.mask=mv.mask|psmask|galmask40
mw.mask=mw.mask|psmask|galmask40
m30.mask=m30.mask|psmask|galmask40
m44.mask=m44.mask|psmask|galmask40
m70.mask=m70.mask|psmask|galmask40
m100.mask=m100.mask|psmask|galmask40
m143.mask=m143.mask|psmask|galmask40
m217.mask=m217.mask|psmask|galmask40
hp.cartview(mw-np.mean(mw[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='WMAP W: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_w.png')
hp.cartview(mv-np.mean(mv[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='WMAP V: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_v.png')
hp.cartview(m30-np.mean(m30[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='Planck 30 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_30.png')
hp.cartview(m44-np.mean(m44[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='Planck 44 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_44.png')
hp.cartview(m70-np.mean(m70[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='Planck 70 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_70.png')
hp.cartview(m100-np.mean(m100[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='Planck 100 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_100.png')
hp.cartview(m143-np.mean(m143[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='Planck 143 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_143.png')
hp.cartview(m217-np.mean(m217[sgppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=True,title='Planck 217 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_217.png')

hp.cartview(hp.smoothing(dvw-np.mean(dvw[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP V-W: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dvw.png')
hp.cartview(hp.smoothing(d3044-np.mean(d3044[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 30-44 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d3044.png')
hp.cartview(hp.smoothing(d4470-np.mean(d4470[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 44-70 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d4470.png')
hp.cartview(hp.smoothing(d70100-np.mean(d70100[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 70-100 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d70100.png')
hp.cartview(hp.smoothing(d100143-np.mean(d100143[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 100-143 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d100143.png')
hp.cartview(hp.smoothing(d143217-np.mean(d143217[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 143-217 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d143217.png')
hp.cartview(hp.smoothing(d44100-np.mean(d44100[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 44-100 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d44100.png')
hp.cartview(hp.smoothing(d70143-np.mean(d70143[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='Planck 70-143 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_d70143.png')

hp.cartview(hp.smoothing(dv44-np.mean(dv44[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP V-Planck 44 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dv44.png')
hp.cartview(hp.smoothing(dv70-np.mean(dv70[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP V-Planck 70 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dv70.png')
hp.cartview(hp.smoothing(dv100-np.mean(dv100[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP V-Planck 100 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dv100.png')
hp.cartview(hp.smoothing(dv143-np.mean(dv143[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP V-Planck 143 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dv143.png')
hp.cartview(hp.smoothing(dw44-np.mean(dw44[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP W-Planck 44 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dw44.png')
hp.cartview(hp.smoothing(dw70-np.mean(dw70[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP W-Planck 70 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dw70.png')
hp.cartview(hp.smoothing(dw100-np.mean(dw100[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP W-Planck 100 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dw100.png')
hp.cartview(hp.smoothing(dw143-np.mean(dw143[sgppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,-90,0],remove_mono=False,title='WMAP W-Planck 143 GHz: SGP $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/sgp_dw143.png')


#NGP

hp.cartview(mw-np.mean(mw[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='WMAP W: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_w.png')
hp.cartview(mv-np.mean(mv[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='WMAP V: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_v.png')
hp.cartview(m30-np.mean(m30[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='Planck 30 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_30.png')
hp.cartview(m44-np.mean(m44[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='Planck 44 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_44.png')
hp.cartview(m70-np.mean(m70[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='Planck 70 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_70.png')
hp.cartview(m100-np.mean(m100[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='Planck 100 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_100.png')
hp.cartview(m143-np.mean(m143[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='Planck 143 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_143.png')
hp.cartview(m217-np.mean(m217[ngppix]),cmap=cmap,min=-250,max=250,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=True,title='Planck 217 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_217.png')

hp.cartview(hp.smoothing(dvw-np.mean(dvw[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP V-W: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dvw.png')
hp.cartview(hp.smoothing(d3044-np.mean(d3044[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 30-44 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d3044.png')
hp.cartview(hp.smoothing(d4470-np.mean(d4470[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 44-70 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d4470.png')
hp.cartview(hp.smoothing(d70100-np.mean(d70100[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 70-100 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d70100.png')
hp.cartview(hp.smoothing(d100143-np.mean(d100143[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 100-143 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d100143.png')
hp.cartview(hp.smoothing(d143217-np.mean(d143217[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 143-217 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d143217.png')
hp.cartview(hp.smoothing(d44100-np.mean(d44100[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 44-100 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d44100.png')
hp.cartview(hp.smoothing(d70143-np.mean(d70143[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='Planck 70-143 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_d70143.png')

hp.cartview(hp.smoothing(dv44-np.mean(dv44[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP V-Planck 44 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dv44.png')
hp.cartview(hp.smoothing(dv70-np.mean(dv70[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP V-Planck 70 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dv70.png')
hp.cartview(hp.smoothing(dv100-np.mean(dv100[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP V-Planck 100 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dv100.png')
hp.cartview(hp.smoothing(dv143-np.mean(dv143[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP V-Planck 143 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dv143.png')
hp.cartview(hp.smoothing(dw44-np.mean(dw44[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP W-Planck 44 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dw44.png')
hp.cartview(hp.smoothing(dw70-np.mean(dw70[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP W-Planck 70 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dw70.png')
hp.cartview(hp.smoothing(dw100-np.mean(dw100[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP W-Planck 100 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dw100.png')
hp.cartview(hp.smoothing(dw143-np.mean(dw143[ngppix]),fwhm=np.radians(.5)),cmap=cmap,min=-25,max=25,lonra=[-20,20],latra=[-20,20],rot=[0,90,0],remove_mono=False,title='WMAP W-Planck 143 GHz: ngp $\pm 10 \deg$',unit='$\mu K$')
savefig('/global/u1/p/peterm/plots/map_comparisons/ngp_dw143.png')





fwhm=2*np.pi/180.
dvwsm2=hp.smoothing(dvw,fwhm=fwhm)
d3044sm2=hp.smoothing(d3044,fwhm=fwhm)
d4470sm2=hp.smoothing(d4470,fwhm=fwhm)
d70100sm2=hp.smoothing(d70100,fwhm=fwhm)
d100143sm2=hp.smoothing(d100143,fwhm=fwhm)
d143217sm2=hp.smoothing(d143217,fwhm=fwhm)
fv44100sm2=hp.smoothing(d44100,fwhm=fwhm)
fv70143sm2=hp.smoothing(d70143,fwhm=fwhm)

dv44sm2=hp.smoothing(dv44,fwhm=fwhm)
dv70sm2=hp.smoothing(dv70,fwhm=fwhm)
dv100sm2=hp.smoothing(dv100,fwhm=fwhm)
dv143sm2=hp.smoothing(dv143,fwhm=fwhm)

dw44sm2=hp.smoothing(dw44,fwhm=fwhm)
dw70sm2=hp.smoothing(dw70,fwhm=fwhm)
dw100sm2=hp.smoothing(dw100,fwhm=fwhm)
dw143sm2=hp.smoothing(dw143,fwhm=fwhm)

fwhm=10*np.pi/180.
dvwsm10=hp.smoothing(dvw,fwhm=fwhm)
d3044sm10=hp.smoothing(d3044,fwhm=fwhm)
d4470sm10=hp.smoothing(d4470,fwhm=fwhm)
d70100sm10=hp.smoothing(d70100,fwhm=fwhm)
d100143sm10=hp.smoothing(d100143,fwhm=fwhm)
d143217sm10=hp.smoothing(d143217,fwhm=fwhm)
fv44100sm10=hp.smoothing(d44100,fwhm=fwhm)
fv70143sm10=hp.smoothing(d70143,fwhm=fwhm)

dv44sm10=hp.smoothing(dv44,fwhm=fwhm)
dv70sm10=hp.smoothing(dv70,fwhm=fwhm)
dv100sm10=hp.smoothing(dv100,fwhm=fwhm)
dv143sm10=hp.smoothing(dv143,fwhm=fwhm)

dw44sm10=hp.smoothing(dw44,fwhm=fwhm)
dw70sm10=hp.smoothing(dw70,fwhm=fwhm)
dw100sm10=hp.smoothing(dw100,fwhm=fwhm)
dw143sm10=hp.smoothing(dw143,fwhm=fwhm)

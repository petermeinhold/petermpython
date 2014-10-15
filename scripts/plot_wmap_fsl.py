sys.path.append('/global/homes/p/peterm/petermpython/paperplots/python/scripts')
sys.path.append('/global/homes/p/peterm/petermpython/paperplots/data')
import cPickle
from setup_matplotlib import *
from matplotlib.colors import ListedColormap
import healpy as hp
cmap = ListedColormap(np.loadtxt("petermpython/paperplots/data/Planck_Parchment_RGB.txt")/255.)


fl=['wmap/wmap_sidelobe_map_W1_9yr_v5.fits',
 'wmap/wmap_sidelobe_map_W2_9yr_v5.fits',
 'wmap/wmap_sidelobe_map_W3_9yr_v5.fits',
 'wmap/wmap_sidelobe_map_W4_9yr_v5.fits',
 'wmap/wmap_sidelobe_map_V1_9yr_v5.fits',
 'wmap/wmap_sidelobe_map_V2_9yr_v5.fits',
 'wmap/wmap_sidelobe_map_V3_9yr_v5.fits']
 
sl_v1=hp.read_map(fl[4],nest=True)
sl_v1r=hp.ud_grade(sl_v1,128,order_in='NEST',order_out='RING')
sl_v1rsm=hp.smoothing(sl_v1r,fwhm=5*np.pi/180.)


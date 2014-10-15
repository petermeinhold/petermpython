import cPickle
from setup_matplotlib import *
from matplotlib.colors import ListedColormap
from glob import glob

cmap = ListedColormap(np.loadtxt("parchment1.dat")/255.)
cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
fl=glob('d*.pkl')
#import newplot


#function to read in a gridded map produced by healpy  (on nersc) along with longitude and latitude
#then produce planck compliant maps using matplotlib instead of healpy
#map should have been saved in a dictionary in a pkl file.
#m = hp.read_map("data/wmap_band_iqumap_r9_7yr_W_v4.fits", 0)

# using directly matplotlib instead of mollview has higher
# quality output, I plan to merge this into healpy

#margins = [0.01, 0.99, 0.99, 0.01]
#hp.mollview(m, min=-1, max=1, unit="mK", title="", xsize=4000, margins=margins)
#hp.graticule()
#plt.savefig("moll.pdf", dpi=200, bbox_inches="tight")


xsize = 2000
ysize = xsize/2.
unit = "$\mu$ K"
vmin = -10; vmax = 10
theta = np.linspace(np.pi, 0, ysize)
phi   = np.linspace(-np.pi, np.pi, xsize)
longitude = np.radians(np.linspace(-180, 180, xsize))
latitude = np.radians(np.linspace(-90, 90, ysize))
PHI, THETA = np.meshgrid(phi, theta)
for f in fl:
    mf=open(f,'rb')
    grid_map=cPickle.load(mf)
    mf.close()
    grid_map=grid_map*1e6
    for width in [8.8]:  #18., 12., 8.8]:
        fig = plt.figure(figsize=(cm2inch(width), cm2inch(width/2.)))
        ax = fig.add_subplot(111,projection='mollweide')

        plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)

        # rasterized makes the map bitmap while the labels remain vectorial
        image = plt.pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=True,cmap=cmap)

        # graticule
        #ax.set_longitude_grid(60)
        #if width < 10:
        #    ax.set_latitude_grid(45)
        #    ax.set_longitude_grid_ends(90)

        # colorbar
        cb = fig.colorbar(image, orientation='horizontal', shrink=.4, pad=0.05, ticks=[vmin, vmax])
        cb.ax.xaxis.set_label_text(unit)
        cb.ax.xaxis.labelpad = -8
        # workaround for issue with viewers, see colorbar docstring
        cb.solids.set_edgecolor("face")

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)


        # remove longitude tick labels
        ax.xaxis.set_ticklabels([])
        # remove horizontal grid
        ax.xaxis.set_ticks([])
        # remove longitude tick labels
        ax.yaxis.set_ticklabels([])
        # remove horizontal grid
        ax.yaxis.set_ticks([])
    
        #plt.grid(True)
        plt.savefig("mollview_"+f[:-4]+"_%dmm.eps" % int(width*10), bbox_inches='tight',dpi=300)
        plt.savefig("mollview_"+f[:-4]+"_%dmm.pdf" % int(width*10), bbox_inches='tight',dpi=300)
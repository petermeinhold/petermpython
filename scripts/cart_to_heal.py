import healpy as hp

def cart_to_heal(cartmap,nside):
    """
    routine to bin a cartesian map into a healpix one
    assume cartmap is a rectangular array
    """
    n_b=cartmap.shape[0]
    n_l=cartmap.shape[1]
    thetas=np.pi-np.arange(n_b)*np.pi/n_b
    phis=np.pi-np.arange(n_l)*np.pi*2./n_l
    outmap=np.zeros(hp.nside2npix(nside),dtype=np.float)
    hitmap=np.zeros(hp.nside2npix(nside))
    for i,phi in enumerate(phis):
        for j,theta in enumerate(thetas):
            iring=hp.ang2pix(nside,theta,phi)
            outmap[iring]=outmap[iring]+cartmap[j,i]
            hitmap[iring]=hitmap[iring]+1
    outmap[hitmap>0]=outmap[hitmap>0]/hitmap[hitmap>0]
    return outmap
    
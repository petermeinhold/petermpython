
import  healpy as hp

def make_gauss_map(nside=64,fwhm=1.0,cpix=0):
    """
    produce healpix map nside=nside with gaussian fwhm=fwhm (degrees) centered on cpix
    """
    npix=hp.nside2npix(nside)
    gmap=np.zeros(npix)
    vcen=hp.pix2vec(nside,cpix)
    sigma=np.radians(fwhm)/(2*np.sqrt(2*log(2.)))
    twosigmasquared=2*sigma**2
    print fwhm/sigma

    for p in arange(npix):
        v=hp.pix2vec(nside,p)
        ang=np.arccos(np.dot(vcen,v))
        gmap[p]=exp(-ang**2/twosigmasquared)
        
    return gmap
    
    
"""set of utility functions for planck analysis"""

import numpy as np
import pyfits as pyfits
import matplotlib.pyplot as plt
from matplotlib import mlab
import time
import os
from scipy import *
from scipy import optimize
from scipy.io import idl
from glob import glob
import pandas as pd
import string
from scipy.io import readsav
import cPickle
import warnings
import sys
sys.path.append('/global/homes/p/peterm/petermpython/paperplots/python/scripts')
sys.path.append('/global/homes/p/peterm/petermpython')
from setup_matplotlib import *
from planckcolors import colombi1_cmap 

warnings.filterwarnings('ignore')
#comment out for local use of planck_util

#import h5py but only if not windows
if os.sys.platform != 'win32':
    import healpy as hp
    import quickring as qr
    import pysqlite2.dbapi2 as sqlite3
    from planck.LFI import LFI
    #import scikits.statsmodels.api as sm
if os.sys.platform=='win32':
    import sqlite3
    
import os
import matplotlib.pyplot as plt

def make_gridmap(inputmap):
    '''
    Function to use correct ct and produce gridmap for paper quality plotting
    '''
    nside = hp.npix2nside(len(inputmap))
    inputmap=hp.ma(inputmap)
    # ratio is always 1/2
    xsize = 2000
    ysize = xsize/2.
    
    theta = np.linspace(np.pi, 0, ysize)
    phi   = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))
    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI)

    # mask
    #if not(inputmap.mask):
    grid_map = inputmap[grid_pix]
    #else:
    #    grid_mask = inputmap.mask[grid_pix]
    #    grid_map = np.ma.MaskedArray(inputmap[grid_pix], grid_mask)
      
    colormaptag = "colombi1_"
    outmap={}
    outmap['cmap'] = colombi1_cmap
    outmap['longitude']=longitude
    outmap['latitude']=latitude
    outmap['data']=grid_map
    return outmap

def save_plot(gridmap,width=8.8,filename='null',figdir='null',grat=True,vmin=-10,vmax=10):
    '''
    function to make and save the paper quality plot, choose width from 18.,12.,8.8
    filename is in PlanckFig_map_<filename>_width.pdf
    '''
    unit = r"$\mathrm{\mu K}$"
    from matplotlib.projections.geo import GeoAxes
    class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
        """Shifts labelling by pi

        Shifts labelling from -180,180 to 0-360"""
        def __call__(self, x, pos=None):
            if x != 0:
                x *= -1
            if x < 0:
                x += 2*np.pi
            return GeoAxes.ThetaFormatter.__call__(self, x, pos)
    
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)))
    # matplotlib is doing the mollveide projection
    ax = fig.add_subplot(111,projection='mollweide')

    # rasterized makes the map bitmap while the labels remain vectorial
    # flip longitude to the astro convention
    image = plt.pcolormesh(gridmap['longitude'][::-1], gridmap['latitude'], gridmap['data'], vmin=vmin, vmax=vmax, rasterized=True, cmap=gridmap['cmap'])

    # graticule
    if grat:
        ax.set_longitude_grid(60)
        ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
        if width < 10:
            ax.set_latitude_grid(45)
            ax.set_longitude_grid_ends(90)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        # remove tick labels
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        # remove grid
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        plt.grid(True)
    # colorbar
    cb = fig.colorbar(image, orientation='horizontal', shrink=.4, pad=0.05, ticks=[vmin, vmax])
    cb.ax.xaxis.set_label_text(unit)
    cb.ax.xaxis.labelpad = -1
    # workaround for issue with viewers, see colorbar docstring
    cb.solids.set_edgecolor("face")

    # remove white space around the image
    # horizontally, vertically the space is removed directly by savefig bbox_inches="tight"
    plt.subplots_adjust(left=0.01, right=0.99)
    plt.savefig(figdir+"PlanckFig_map_" + filename + "_%dmm.pdf" % int(width*10), bbox_inches='tight', pad_inches=0.02)
    
def savef(path, ext='png', close=True, verbose=False):
    """
    Save a figure from pyplot. (from http://www.jesshamrick.com/2012/09/03/saving-figures-from-pyplot/)
    Parameters
    ----------
    path : string
    The path (and filename, without the extension) to save the
    figure to.
    ext : string (default='png')
    The file extension. This must be supported by the active
    matplotlib backend (see matplotlib.backends module). Most
    backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.  
    close : boolean (default=True)
    Whether to close the figure after saving. If you want to save
    the figure multiple times (e.g., to multiple formats), you
    should NOT close it in between saves or you will have to
    re-plot it.
    verbose : boolean (default=True)
    Whether to print information about when and where the image
    has been saved.
    """
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
    # The final path to save to
    savepath = os.path.join(directory, filename)
    if verbose:
        print("Saving figure to '%s'..." % savepath),
    # Actually save the figure
    plt.savefig(savepath)
    # Close it
    if close:
        plt.close()
    if verbose:
        print("Done")
    
def get_lfi_dx11_mask(nside):
    maskfile=open('/global/homes/p/peterm/masks/dx11_lfi_total_mask_2048.pkl','rb')
    tmask=cPickle.load(maskfile)
    maskfile.close()
    tmask=degrade_mask(tmask,nside_out=nside)
    return tmask
    
def get_ffp8_cls(freq,s1,s2,hr1,hr2):
    #find and read the pickled CL list from FFP8 noise sims
    topdir='/global/homes/p/peterm/ffp8_noise_cls/'
    pklfilename=topdir+'ffp8_noise_null_cls_'+freq+s1+hr1+s2+hr2+'.pkl'
    pklfile=open(pklfilename,'rb')
    cls=cPickle.load(pklfile)
    pklfile.close()
    return cls
    
def get_ffp8_mean_fit_cls(freq,q1,q2,s1,s2,hr1,hr2,det=''):
    #find and read the pickled CL list from FFP8 noise sims
    
    topdir='/global/homes/p/peterm/ffp8_noise_cls/'
    pklfilename=topdir+'cl_fit_ffp8_noise_null_cls_'+freq+det+q1+s1+hr1+q2+s2+hr2+'.pkl'
    pklfile=open(pklfilename,'rb')
    cls=cPickle.load(pklfile)
    pklfile.close()
    return cls
    
def fitffpcls(cls):
    #function takes list of TT,EE,BB cls from ffp sim nulls, histograms and fits asymmetric gaussian errors
    cls=np.array(cls)
    ncls=len(cls[0,0,:])
    allfitcls={}
    for p,plabel in zip([0,1,2],['TT','EE','BB']):
        fitcls={}
        meancls=np.zeros(ncls,dtype=float)
        lowsigcls=np.zeros(ncls,dtype=float)
        hisigcls=np.zeros(ncls,dtype=float)
        for lbin in  range(ncls):
            h=histogram(cls[:,p,lbin],bins=20)
            st=np.std(cls[:,p,lbin])
            mn=np.mean(cls[:,p,lbin])
            pstart=[len(cls[:,p,lbin])/10.,mn,st,st]
            test=fit_asymmetric_gaussian(h[1][1:],h[0],pstart=pstart)
            meancls[lbin]=test['center']
            lowsigcls[lbin]=test['sigmalow']
            hisigcls[lbin]=test['sigmahigh']
        fitcls['l']=arange(1,ncls+1)
        fitcls['mean']=meancls
        fitcls['lowerr']=lowsigcls
        fitcls['higherr']=hisigcls
        allfitcls[plabel]=fitcls

    return allfitcls

def read_and_diff_files(f1,f2,nside=None,tmask=None,corr1=None,corr2=None,return_map=False):
    colnames=['I_Stokes','Q_Stokes' ,'U_Stokes' ,'Hits    ' ,'II_cov  ' ,'IQ_cov  ' ,'IU_cov  ' ,'QQ_cov  ' ,'QU_cov  ' ,'UU_cov  ' ]
       
    m1={}
    m2={}
    for i in range(3):
        mm1=hp.ma(hp.read_map(f1,[i],verbose=False))
        m1[colnames[i].strip()]=mm1
        mm2=hp.ma(hp.read_map(f2,[i],verbose=False))
        m2[colnames[i].strip()]=mm2

    keys=m1.keys()
    keys
    # Generate differences
    mdiff={}
    for key in ['I_Stokes','Q_Stokes','U_Stokes']:
        if tmask is None:
            tmask= m1[key].mask | m2[key].mask
        tmask=tmask | m1[key].mask | m2[key].mask
        
    if nside:
        tmask=degrade_mask(tmask,nside_out=nside)        
    for key in ['I_Stokes','Q_Stokes','U_Stokes']:
        if corr1 is None:
            mdiff[key]=hp.ma(m1[key]-m2[key])
        if corr1:
            nside_map=hp.npix2nside(len(m1[key]))
            corr1udg=hp.ud_grade(corr1[key],nside_map)
            corr2udg=hp.ud_grade(corr2[key],nside_map)
            mdiff[key]=hp.ma((m1[key]-corr1udg)-(m2[key]-corr2udg))
    skyfrac=1-float(tmask.sum())/len(tmask)            
    # Get MC noise per pixel IQU
    mdiffd=[]
    for key in ['I_Stokes','Q_Stokes','U_Stokes']:
        if nside==None:
            mkey=mdiff[key]
            mkey.mask=tmask
        if nside!=None:
            mkey=hp.ud_grade(mdiff[key],nside_out=nside)    
            mkey.mask=tmask
        mdiffd.append(mkey)
    cldata=hp.anafast(mdiffd)
    cldata_out=[]
    for cl in cldata:
        cldata_out.append(cl/skyfrac)
    cldatad={}
    for i,spec in enumerate(['TT','EE','BB','TE','TB','EB']):
        cldatad[spec]=cldata_out[i]
    print 'skyfrac ',skyfrac
    if return_map is False:
        return cldatad
    if return_map is True:
        return cldatad,mdiffd

def degrade_mask(inmask,nside_out=256):
    #stupid function to degrade a mask by making a map and degrading that
    m=hp.ma(np.ones(len(inmask)))
    m.mask=inmask
    mdg=hp.ud_grade(m,nside_out=nside_out)
    return mdg.mask
    
    
def asgaussian(x,p):
    """
    define an asymmetric gaussian for error fitting of slices
    parameters p(4): (amplitude,center, sigmaleft,sigmaright)
    """
    ind1=x<=p[1]
    ind2=x>p[1]
    x1=x[ind1]
    x2=x[ind2]
    yleft=p[0]*np.exp(-(p[1]-x1)**2/(2*p[2]**2))
    yright=p[0]*np.exp(-(p[1]-x2)**2/(2*p[3]**2))
    y=np.hstack((yleft,yright))
    return y

def asgresiduals(p,y,x):
    err=y-asgaussian(x,p)
    return err
    
def fit_asymmetric_gaussian(x,y,pstart=None):
    """
    function to fit asymmetric gaussian
    """
    xmax=x[y==np.max(y)][0]
    if pstart==None:
        pstart=[np.max(y),xmax,1,1]
        if np.max(y)==y[0]:
            pstart=[y[0],x[0]-1,0,1]
        if np.max(y)==y[-1]:
            pstart=[y[-1],x[-1]+1,1,0]
    pbest=optimize.leastsq(asgresiduals,pstart,args=(y,x))
    params={}
    params['amplitude']=pbest[0][0]
    params['center']=pbest[0][1]
    params['sigmalow']=np.abs(pbest[0][2])
    params['sigmahigh']=np.abs(pbest[0][3])
    params['fit']=asgaussian(x,pbest[0])
    return params
    
def read_like_slices(prefix='slices_pm',xspec='TE',directory='.'):
    """
    function to read comm_like_tools likelihood slices, fit ell bins and return a dictionary with bins, means, asymmetric errors
    prefix = slicefileprefix to choose slice files
    xspec chooses which to read in (must ha
    ve been included in bins.dat)
    directory defaults to current directory
    """
    fl=glob(directory+'/'+prefix + '_'+ xspec+'*.dat')
    fl.sort()
    if len(fl)>0:
        spec={}
        ellbins=[]
        values=[]
        lowval=[]
        hival=[]
        for f in fl:
            lowel=np.int(f[-13:-9])
            hiel=np.int(f[-7:-4])
            d=genfromtxt(f)
            if len(d[:,0]) > 3:
                fit=fit_asymmetric_gaussian(d[:,0],d[:,1])
                values.append(fit['center'])
                sigmalow=fit['sigmalow']
                sigmahi=fit['sigmahigh']
                #If fit failed on one side, use a symmetric eror bar
                #if sigmahi > 5*sigmalow:
                #    sigmahi=sigmalow
                #if sigmalow > 5*sigmahi:
                #    sigmalow=sigmahi
                lowval.append(sigmalow)
                hival.append(sigmahi)
                ellbins.append(np.int(np.mean([lowel,hiel])))
        spec['ellbins']=np.array(ellbins)
        spec['values']=np.array(values)
        spec['lowval']=np.array(lowval)
        spec['hival']=np.array(hival)
    else:
        print 'No matching files'
        spec=None
    return spec
    

def od2pid_range(od,ucds=None,key='LFI27M'):
    """
    function returns 2 element range min and max pID of the input OD. If UCDS is None, reads it in
    """
    if ucds==None:
        ucds=get_ucds(freq=30)
    pidrange=[np.min(ucds[key]['pID'][ucds[key]['od_int']==od]),np.max(ucds[key]['pID'][ucds[key]['od_int']==od])]
    return pidrange

def logbin(freq,ps,nbins):
    """
    function to bin in logarithmic frequency intervals
    output freq and ps arrays will be nbins long
    """
    logfreqs=logspace(np.log2(np.min(freq[freq>0])),np.log2(np.max(freq)),num=nbins,endpoint=True,base=2)
    
    outps=[]
    outfreqs=[]
    lastfreq=0.
    for logfreq in logfreqs:
        outps.append(np.mean(ps[np.logical_and(freq > lastfreq,freq<=logfreq)]))
        #outfreqs.append(np.mean(freq[np.logical_and(freq >= lastfreq,freq<=logfreq)]))
        outfreqs.append((logfreq+lastfreq)/2.)
        lastfreq=logfreq
    return np.array(outfreqs),np.array(outps)
    
def get_mapname(caltype='DCRS',channel=None,quad=[20,21],survey=None,subchunk=None,pol=0):
    #quick util to read map from Andreas folders
    #parse to name:
    basedir='/global/homes/z/zonca/m/maps/out/'
    if quad!=None:
        chname=''
        for horn in quad:
            for rad in ['M','S']:
                chname=chname+'LFI'+str(horn)+rad+'_'
    if quad==None:
        if channel== None:
            chname='LFI18M_'
        if channel!=None:
            chname=channel+'_'
    if survey==None:
        surveyname='nominal'
    if survey!=None:
        surveyname='survey'+str(survey)
    if subchunk==None:
        chunk=''
    if subchunk!=None:
        chunk='_subchunk_'+str(subchunk)+'_of_2'
    fname=basedir+'map_'+caltype+'_'+chname+surveyname+'_'+'1min'+chunk+'.fits'
    ffname=glob(fname)[-1]
    #m=hp.ma(hp.read_map(ffname),pol)
    return ffname
    
    
def get_rel_scales(psmap,reverse=False,usevar=True):
    """
    function to find linear fit params among adjacent
    phase binned rings in a pseudomap
    returns array of slopes and offsets
    row_i ~ fitpar1*row_(i+1) + fitpar2
    """
    nrows=np.shape(psmap)[0]
    fitpars1=[]
    fitpars2=[]
    for i in range(nrows-1):
        if usevar==True:
            fp1=np.var(psmap[i,:])/np.var(psmap[i+1,:])
            fp2=0
        if usevar==False:
            if reverse==True:
                fp1,fp2=polyfit(psmap[i,:],psmap[i+1,:],1)
            if reverse==False:
                fp1,fp2=polyfit(psmap[i+1,:],psmap[i,:],1)
        fitpars1.append(fp1)
        fitpars2.append(fp2)
    fitpars1=np.array(fitpars1)
    fitpars2=np.array(fitpars2)
    return fitpars1,fitpars2


def maskmap(m):
    m_out = hp.ma(m)
    m_out.mask = m == 0
    return m_out.filled()
    
def mapsum (inputfilelist,freq=70,wmapmask=True,psmaskmask=True,single=False):
    """function to read in a list of maps and sum them using
    a union mask. Default options to include wmap mask and
    pt source mask. Assume IQU inputs, read separately, mask
    sum and into list default freq(for psmask) =70
    if single==true, assume input maps are T only, and second column is Hits, use this to fill all the map
    """
    #read one map, just to find out size
    m=hp.ma(hp.read_map(inputfilelist[0]),0)
    nside=hp.npix2nside(len(m))
    tmask=m.mask
    tsum=np.zeros(len(m))
    hsum=np.zeros(len(m))

    if wmapmask==True:
        wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
        tmask=tmask|wmask
    if psmaskmask==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        tmask=tmask|psmask
    nmaps=len(inputfilelist)
    if single==False:
        qsum=np.zeros(len(m))
        usum=np.zeros(len(m))
        for f in inputfilelist:
            t=hp.ma(hp.read_map(f,0))
            q=hp.ma(hp.read_map(f,1))
            u=hp.ma(hp.read_map(f,2))
            tmask=tmask|t.mask|q.mask|u.mask
            t.mask=tmask
            q.mask=tmask
            u.mask=tmask
            tsum=tsum+t/nmaps
            qsum=qsum+q/nmaps
            usum=usum+u/nmaps
        outmap=[tsum,qsum,usum]
    if single==True:
        for f in inputfilelist:
            t=hp.ma(hp.read_map(f,0))
            h=hp.ma(hp.read_map(f,1))
            tmask=tmask & t.mask
            t.mask=tmask
            hsum=hsum+h
            tsum=tsum+t*h
        tsum=tsum/hsum
        outmap=[tsum,hsum]
    return outmap

def make_bmask(nside=1024,bmax=30):
    """
    function to provide simple galactic latitude cut mask
    returns masked array of ones and zeros where mask is
    """
    npix=hp.nside2npix(nside)
    blk=hp.ma(np.arange(npix))
    theta,phi=hp.pix2ang(nside,blk)
    b=180*(theta-np.pi/2.)/np.pi
    blk=hp.ma(np.ones(npix))
    blk.mask=np.abs(b)<bmax
    blk[blk.mask]=0
    blk=hp.ma(maskmap(blk))
    return blk
    
def get_ring_for_spin_vec( vec, npts=500 ):
    """ returns a [0:3,0:nphi] length array containing nphi
        points along the ring which intersects the unit sphere
        and the plane which is normal to vec."""
    
    # get a good vector which is not colinear-specialized for galactic coord
    uvec = np.array([-0.09647625,  0.86228144,  0.49715496])

    # get two basis vectors which span the ring.
    bas1  = np.cross( vec, uvec )
    bas1 *= np.sqrt( (1.0 - np.dot(vec, vec)) / np.dot(bas1, bas1) )

    bas2  = np.cross(vec, bas1)
    bas2 *= np.sqrt( (1.0 - np.dot(vec, vec)) / np.dot(bas2, bas2) )

    # form a curve for the ring.
    ang = np.linspace(0, 2.*np.pi, npts)
    xs  = vec[0] + bas1[0]*np.cos(ang) + bas2[0]*np.sin(ang)
    ys  = vec[1] + bas1[1]*np.cos(ang) + bas2[1]*np.sin(ang)
    zs  = vec[2] + bas1[2]*np.cos(ang) + bas2[2]*np.sin(ang)

    return xs, ys, zs

def rescan_to_psmap(inputmap,startring=0,stopring=5500,stepring=10,nphi=500,ang=85):
    #function to make pseudomap from (assumed galactic) map
    x=arange(nphi)
    nside=hp.npix2nside(len(inputmap))
    vecs=qr.get_hfi_ring_spin_axes()* np.cos( ang * np.pi / 180. )
    psmap=[]
    inputmap=hp.ma(inputmap)
    inputmap.fill_value=inputmap.mean()
    inmap=inputmap.filled()
    for ring in range(startring,stopring,stepring):
        vec = vecs[ring] 
        xs, ys, zs = get_ring_for_spin_vec(vec, npts=nphi)
        tht   = np.pi*0.5 - np.arcsin(zs) #latitude
        phi   = np.arctan2(ys, xs)        #longitude
        tpix=hp.ang2pix(nside,tht,phi)
        psmap.append(inmap[tpix])
    psmap=np.array(psmap)
    return psmap
    
def import_csv_cals(infile):
    """
    function to read in ascii file, assume first line has strings naming columns
    """
    f=open(infile,'r')
    b=f.readlines()
    f.close()
    names=np.str.split(b[0],' ')
    ncol=len(names)
    pid=[]
    gainm=[]
    gains=[]
    for row in b[1:]:
        srow=np.str.strip(row,'\n')
        rrow=np.str.split(srow,' ')
        pid.append(np.float(rrow[0]))
        gainm.append(np.float(rrow[1]))
        gains.append(np.float(rrow[2]))
    pid=np.array(pid)
    gainm=np.array(gainm)
    gains=np.array(gains)
    outgains=np.vstack([pid,gainm,gains])
    return outgains
        
        

def read_bandpass_weights(filename='../frequency_bandpasses/bandpassweights_2.csv'):
    """
    function to import and reformat 30 GHz bandpass weights from csv file to flat array
    """
    f=open(filename,'r')
    b=f.readlines()
    f.close()
    bb=np.zeros([7,5],dtype=float)
    for i in range(7):
        bb[i,:]=np.array(np.str.split(b[i+1],','))
    
    bbb=bb[:,1:].T
    bbb=bbb.reshape(28)
    return bbb


def plot_fsl_model_resid(datadiffmap,fsldiffmaplist,mask,weights,norm=1.,min=-20e-6,max=5e-6,title='Residual'):
    """
    function to mollview plot datadiffmap-model where model is
    calculated by summing all fsl sims with 28 row weights vector
    norm is applied after normalizing weight vector to 1 total.
    """
    weights=norm*weights/np.sum(weights)
    model=hp.ma(np.sum(fsldiffmaplist*weights,axis=1))
    model.mask=mask
    pl1=hp.mollview(model,min=min,max=max,title='Model')
    plt.show()
    rmsdata=np.std(datadiffmap[mask==False])
    rmsres=np.std(datadiffmap[mask==False]-model[mask==False])
    pl2=hp.mollview(datadiffmap-model,min=min,max=max,title=title+'RMS data/resid '+np.str(rmsdata)+' / '+np.str(rmsres))
    plt.show()
    return
    
def get_surv_diff_chan(chan1='27M',surv1=1,chan2='27M',surv2=2,mask_ps=True,nside=1024,fwhm=2.0):
    """
    function to retrieve two single channel survey maps, common mask, and difference then smooth 
    chan is of form '27M', maps are fits, fwhm is in degrees, if zero or negative don't smooth at all, default 2 degrees
    """
    
    pol='I'
    freq=70
    if np.int(chan1[0:2])>23:
        freq=44
    if np.int(chan1[0:2])>26:
        freq=30

    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    #base='/global/u1/p/peterm/dx9_single_horn_survey/'
    base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/channels/'
    mapname='*'+np.str(chan1)+'*survey_'+str(surv1)+'.fits'
    fl=glob(base+mapname)
    m1=hp.ma(hp.ud_grade(maskmap(hp.read_map(fl[0])),nside)) 
        
    mapname='*'+np.str(chan2)+'*survey_'+str(surv2)+'.fits'
    fl=glob(base+mapname)
    m2=hp.ma(hp.ud_grade(maskmap(hp.read_map(fl[0])),nside)) 
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=totalmask
    if fwhm<=0:
        difsm=dif
    return difsm

def get_surv_diff_fname(file1,file2,pol='I',nside=1024,freq=30,mask_ps=True,fwhm=2.0):
    """
    function to retrieve two survey maps, common mask, and difference then smooth to 2 degrees
    input the filenames themselves need freq to get correct psmask
    fwhm is in degrees, if zero or negative don't smooth at all, default 2 degrees
    """
    npol=0
    if pol=='Q':
        npol=1
    if pol=='U':
        npol=2
    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    m1=hp.ma(hp.ud_grade(maskmap(hp.read_map(file1,npol)),nside)) 
        
    m2=hp.ma(hp.ud_grade(maskmap(hp.read_map(file2,npol)),nside)) 
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=totalmask
    if fwhm<=0:
        difsm=dif
    return difsm
    
def surv_diff(surveydict,ss1=1,ss2=2,nside=128,freq=70,mask_ps=True,fwhm=10.0):
    """
    function to make differences among 5 surveys from andreas interactive destriper/binner. uses common masks smooth to 10 degrees
    fwhm is in degrees, if zero or negative don't smooth at all, default 10 degrees
    """
    m1=hp.ma(np.array(surveydict[ss1]))
    m2=hp.ma(np.array(surveydict[ss2])) 
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside,order_out='NEST')))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    dif.mask |= np.isnan(dif)
    dif=hp.ud_grade(dif,nside,order_in='NEST',order_out='RING')
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=totalmask
    if fwhm<=0:
        difsm=dif
    return difsm
def get_surv_diff(chan1=30,surv1=1,cal1='dx9',pol1='I',chan2=30,surv2=2,cal2='dx9',pol2='I',mask_ps=True,nside=1024,fwhm=2.0):
    """
    function to retrieve two survey maps, common mask, and difference then smooth to 2 degrees
    if using ouptut of DST, cal is of form 'dx8s_30', denotes frequency , maps are H5, if using
    dpc maps, cal is of form 'dx9', chan is of form '30', maps are fits
    fwhm is in degrees, if zero or negative don't smooth at all, default 2 degrees
    """
    freq=chan1
    npol1=0
    if pol1=='Q':
        npol1=1
    if pol1=='U':
        npol1=2
    npol2=0
    if pol2=='Q':
        npol2=1
    if pol2=='U':
        npol2=2

    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    if ((cal1 != 'dx9') & (cal1 != 'dx8') & (cal1 != 'dx7')):
        base='/global/scratch/sd/planck/user/zonca/resdst/'
        nside=1024
        dirname=cal1
        mapname='/surv'+np.str(surv1)+'/map_ch'+np.str(chan1)+'.h5'
        m1= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+mapname,mode='r')['data'][pol1][:]),nside,order_in='NESTED',order_out='RING'))
    if ((cal1 == 'dx9')|(cal1=='dx8')|(cal1=='dx7')):
        base='/global/project/projectdirs/planck/data/mission/DPC_maps/'+cal1+'/lfi/'
        if surv1!=0:
            mapname='LFI_'+np.str(chan1)+'*survey_'+str(surv1)+'.fits'
        if surv1==0:
            mapname='LFI_'+np.str(chan1)+'*_full.fits'
        pol=pol1
        fl=glob(base+mapname)
        m1=hp.ma(hp.ud_grade(maskmap(hp.read_map(fl[0],npol1)),nside)) 
        
    if ((cal2 != 'dx9') & (cal2 != 'dx8') & (cal2 != 'dx7')):
        base='/global/scratch/sd/planck/user/zonca/resdst/'
        nside=1024
        dirname=cal2
        mapname='/surv'+np.str(surv2)+'/map_ch'+np.str(chan2)+'.h5'
        m2= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+mapname,mode='r')['data'][pol2][:]),nside,order_in='NESTED',order_out='RING'))
    if ((cal2 == 'dx9')|(cal2=='dx8')|(cal2=='dx7')):
        base='/global/project/projectdirs/planck/data/mission/DPC_maps/'+cal1+'/lfi/'
        if surv2!=0:
            mapname='LFI_'+np.str(chan2)+'*survey_'+np.str(surv2)+'.fits'
        if surv2==0:
            mapname='LFI_'+np.str(chan2)+'*_full.fits'
        pol=pol2
        fl=glob(base+mapname)
        m2=hp.ma(hp.ud_grade(maskmap(hp.read_map(fl[0],npol2)),nside))
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    dif.mask |= np.isnan(dif)
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=totalmask
    if fwhm<=0:
        difsm=dif
    return difsm

def get_surv_diff_fsl(chan=27,pol=0,rad='M',mask_ps=True,nside=512,maptype='destriped',subband=1,dx9template=False,dipole=False,fwhm=2.0):
    """
    function to retrieve fsl simulation survey maps, common mask, and difference(ss1-ss2) then smooth to 2 degrees
    hardwired to location in peterm home on carver .  pol= integer: 0='I',1='Q',2='U',subband=int from 0 to 6
    fwhm is in degrees, if zero or negative don't smooth at all, default 2 degrees
    """
    freq=70
    if chan>23:
        freq=44
    if chan>26:
        freq=30
    polstring=''
    if pol!= 0:
        polstring='_pol'
    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    if dx9template==False:
        if dipole==False:
            if maptype=='binned':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss1_512_binned.fits'
            if maptype=='destriped':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss1_512.fits'
        if dipole ==True:
            if maptype=='binned':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/dipole/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss1_512_binned_dipole.fits'
            if maptype=='destriped':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/dipole/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss1_512_dipole.fits'
            
    if dx9template==True:
        if maptype=='binned':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
            mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss1_512_binned_dx9.fits'
        if maptype=='destriped':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
            mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss1_512_dx9.fits'

    m1=hp.ma(maskmap(hp.ud_grade(hp.read_map(mapname,pol),nside))) 

    if dx9template==False:
        if dipole==False:
            if maptype=='binned':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss2_512_binned.fits'
            if maptype=='destriped':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss2_512.fits'
        if dipole ==True:
            if maptype=='binned':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/dipole/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss2_512_binned_dipole.fits'
            if maptype=='destriped':
                base='/global/u1/p/peterm/sims_30ghz_bandpass/dipole/'
                mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss2_512_dipole.fits'

    if dx9template==True:
        if maptype=='binned':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
            mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss2_512_binned_dx9.fits'
        if maptype=='destriped':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
            mapname=base+'madam_slb_LFI'+np.str(chan)+rad+polstring+'_f'+np.str(subband)+'_ss2_512_dx9.fits'

    m2=hp.ma(maskmap(hp.ud_grade(hp.read_map(mapname,pol),nside))) 
    
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=totalmask
    if fwhm<=0:
        difsm=dif
    return difsm
    
def get_surv_diff_fsl_pol(freq=30,pol=0,mask_ps=True,nside=512,maptype='destriped',subband=1,dx9template=True,dipole=False,fwhm=2.0):
    """
    function to retrieve fsl simulation survey maps, common mask, and difference(ss1-ss2) then smooth to 2 degrees
    hardwired to location in peterm home on carver .  pol= integer: 0='I',1='Q',2='U',subband=int from 0 to 6
    madam_slb_30GHz_pol_f1_ss1_dx9_512.fits
    fwhm is in degrees, if zero or negative don't smooth at all, default 2 degrees
    """
    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    if dx9template==False:
        if maptype=='binned':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss1_512_binned.fits'
        if maptype=='destriped':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss1_512.fits'
            
    if dx9template==True:
        if maptype=='binned':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss1_dx9_512_binned.fits'
        if maptype=='destriped':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss1_dx9_512.fits'

    m1=hp.ma(maskmap(hp.ud_grade(hp.read_map(mapname,pol),nside))) 

    if dx9template==False:
        if maptype=='binned':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss2_512_binned.fits'
        if maptype=='destriped':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss2_512.fits'
            
    if dx9template==True:
        if maptype=='binned':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/binned/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss2_dx9_512_binned.fits'
        if maptype=='destriped':
            base='/global/u1/p/peterm/sims_30ghz_bandpass/destriped/'
            mapname=base+'madam_slb_'+np.str(freq)+'GHz_pol_f'+np.str(subband)+'_ss2_dx9_512.fits'

    m2=hp.ma(maskmap(hp.ud_grade(hp.read_map(mapname,pol),nside))) 
    
    totalmask=m1.mask|m2.mask
    if mask_ps==True:
        #psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))
        totalmask=m1.mask | m2.mask | psmask
    dif=(m1-m2)/2.
    dif.mask=totalmask
    if fwhm>0:
        difsm=hp.ma(hp.smoothing(dif.filled(),fwhm*np.pi/180.))
        difsm.mask=totalmask
    if fwhm<=0:
        difsm=dif
    return difsm


def get_dst_map(chan=None,surv='surv1',cal='dx8s_30',pol='I'):
    """
    function to retrieve a survey map
    """
    nside=1024
    dirname=cal
    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    base='/global/scratch/sd/planck/user/zonca/resdst/'
    if chan==None:
        mapname='/'+surv+'/map.h5'
        m= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+mapname,mode='r')['data'][pol][:]),nside,order_in='NESTED',order_out='RING'))    
    if chan!=None:
        mapname='/'+surv+np.str(surv1)+'/map_ch'+np.str(chan1)+'.h5'
        m= hp.ma(hp.ud_grade(maskmap(h5py.File(base+dirname+mapname,mode='r')['data'][:]),nside,order_in='NESTED',order_out='RING')) 
    return m
    

def get_dx10_map(chan='070',surv='yr1',half='',mask_ps=True,mask_gal=True,nside=1024,pol=None,eb=False,fwhm=None,dx11_pre=False):
    """
    function to retrieve survey map from DX10 releases. chan can be '070-19M' or '070'
    survey can be 'yr1', 'full','survey_1' ...., half can be none or 1 or 2
    """
    freq=np.int(chan[:3])
    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    if np.float(freq)<90:
        base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/lfi/'
        file_nside=1024
        lfi=True
        hfi=False
    if np.float(freq)>90:
        if dx11_pre==False:
            base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/hfi/official/'
        if dx11_pre==True:
            base='/project/projectdirs/planck/data/mission/DPC_maps/dx11_pre/hfi/extra/standard_dpc_gains/FREQ/'
        file_nside=2048
        hfi=True
        lfi=False

    if half==2: 
        half='_ringhalf_2'
    if half==1: 
        half='_ringhalf_1'
    
    mapname='*SkyMap*'+chan+'_'+str(file_nside)+'*'+surv+half+'.fits'
    if pol=='I':
        poln=0
    if pol=='Q':
        poln=1
    if pol=='U':
        poln=2
    

    print base+mapname
    fl=glob(base+mapname)
    print fl
    if pol == None:
        m= hp.ma(hp.ud_grade(hp.read_map(fl[0],field=(0,1,2)),nside))
        totalmask= m[1] == 0
    if pol != None:
        m= hp.ma(hp.ud_grade(hp.read_map(fl[0],field=poln),nside))
        totalmask= m == 0
    if mask_ps==True:
        if hfi:
            psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map('/project/projectdirs/planck/data/releases/HFI_PR1/HFI_Mask_PointSrc_2048_R1.10.fits'), nside)))
        if lfi:
            psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/releases/LFI_PR1/' + 'LFI_MASK_0%d-ps_*.fits' % freq)[0]), nside)))

        totalmask=totalmask | psmask
    if mask_gal==True:
        galmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map('/project/projectdirs/planck/data/releases/LFI_PR1/COM_MASK_gal-06_2048_R1.00.fits'), nside)))
        totalmask=totalmask | galmask
    if pol == None:
        for map in m:
            map.mask = totalmask        
    if eb:
        #convert tqu to teb using alm
        nside=hp.npix2nside(len(m[0]))
        alm=hp.map2alm(m,lmax=nside*4)
        m=[]
        m.append(hp.ma(hp.alm2map(alm[0],nside,lmax=nside*4)))
        m.append(hp.ma(hp.alm2map(alm[1],nside,lmax=nside*4)))
        m.append(hp.ma(hp.alm2map(alm[2],nside,lmax=nside*4)))
    if pol==None:
        for map in m:
            map.mask=totalmask
    if fwhm != None:
        mcopy=copy(m)
        m=[]
        for map in mcopy:
            mapsm=hp.ma(hp.smoothing(map,fwhm=fwhm))
            mapsm.mask=totalmask
            m.append(mapsm)
    return m


def get_dx_map(chan='70',surv='yr1_ringhalf_1',cal='dx9',mask_ps=False,pol='I'):
    """
    function to retrieve survey map from DX releases
    """
    
    freq=np.int(chan)
        
    def maskmap(m):
        m_out = hp.ma(m)
        m_out.mask = m == 0
        return m_out.filled()
    base='/global/project/projectdirs/planck/data/mission/DPC_maps/'
    dirname=cal.lower()
    calname=cal.upper()
    nside=1024
    mapname='/lfi/'+'LFI_*'+chan+'*'+surv+'.fits'

    if pol=='I':
        poln=0
    if pol=='Q':
        poln=1
    if pol=='U':
        poln=2
    print base+dirname+mapname
    fl=glob(base+dirname+mapname)
    print fl
    m= hp.ma(maskmap(hp.read_map(fl[0])))
    if mask_ps==True:
        psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
        totalmask=m.mask | psmask
    m.mask=totalmask
    return m
    
    

#def get_ucds(db="/project/projectdirs/planck/user/zonca/issues/deltav/ucds-dx9-delta.db",freq=30):
def get_ucds(db="/global/homes/p/peterm/ucds-dx10.db",freq=30):
    MS = 'MS'
    #conn = sqlite3.connect("/project/projectdirs/planck/user/zonca/issues/deltav/ucds-dx9.db")
    conn = sqlite3.connect(db)
    c= conn.cursor()
    ch_data_dtype = np.dtype(
        [('pID',np.long), ('obt_start',np.long),('obt_end',np.long),('od_int',np.int),('vref',np.double),('vsky',np.double),('vref0',np.double),('vref1',np.double), 
        ('vsky0',np.double),('vsky1',np.double), 
        ('gain',np.double), ('nominal_gain',np.double), ('dipoleT',np.double), 
        ('l14k',np.double), ('ts1l',np.double), ('ts2l',np.double), ('ts3l',np.double), 
        ('ts4l',np.double), ('ts5l',np.double), ('ts6l',np.double), ('ts1r',np.double), 
        ('ts2r',np.double), ('ts3r',np.double), ('ts4r',np.double), ('ts5r',np.double), 
        ('ts6r',np.double), ('beusvc',np.double), ('beusci1',np.double), ('beusci2',np.double), 
        ('beusci3',np.double), ('beusci4',np.double), ('lbem1',np.double), ('lbem2',np.double), 
        ('lfem1',np.double), ('lfem2',np.double), ('ldaq1',np.double), ('rbem1',np.double), 
        ('rbem2',np.double), ('rfem1',np.double), ('rfem2',np.double), ('rdaq1',np.double)])
    f = LFI().f[freq]
    #construct set of hskp data vectors to fit to:
    all_data = {}
    for ch in f.ch:
        horn = ch.horn
        rad = ch.n
        chtag = "LFI%d%s" % (horn, MS[rad])
        c.execute(' '.join([ "select sci.pointingID as pID,",
             "pointing_info.start_time as obt_start,",
             "pointing_info.end_time as obt_end,",
             "pointing_info.od_int as od_int,",
             "sci.mean_ref as vref,",
             "sci.mean_sky as vsky,",
             "sci_d0.mean_ref as vref0,",
             "sci_d1.mean_ref as vref1,",
             "sci_d0.mean_sky as vsky0,",
             "sci_d1.mean_sky as vsky1,",
             "rg.gain as gain,",
             "sg.gain as nominal_gain,",
             "dipole.deltaT as dipoleT,",
             "hfihk.l1_a5 as l14k,",
             "lfihk.feu_left_side_wall as ts1l,",
             "lfihk.feu_cold_plate_left_inner as ts2l,",
             "lfihk.feu_cold_plate_right_inner as ts3l,",
             "lfihk.feu_left_bottom_fh25 as ts4l,",
             "lfihk.feu_cold_plate_far_left as ts5l,",
             "lfihk.feu_cone_left_part as ts6l,",
             "lfihk.feu_right_bottom_fh28 as ts1r,",
             "lfihk.feu_cone_right_part as ts2r,",
             "lfihk.feu_right_side_wall as ts3r,",
             "lfihk.feu_cold_plate_far_right as ts4r,",
             "lfihk.feu_fh28_flange as ts5r,",
             "lfihk.feu_right_bottom_fh26 as ts6r,",
             "lfihk.beu_service as beusvc,",
             "lfihk.beu_science1 as beusci1,",
             "lfihk.beu_science2 as beusci2,",
             "lfihk.beu_science3 as beusci3,",
             "lfihk.beu_science4 as beusci4,",
             "lfihk.beu_l_bem1 as lbem1,",
             "lfihk.beu_l_bem2 as lbem2,",
             "lfihk.beu_l_fem1 as lfem1,",
             "lfihk.beu_l_fem2 as lfem2,",
             "lfihk.beu_r_bem1 as rbem1,",
             "lfihk.beu_r_bem2 as rbem2,",
             "lfihk.beu_r_fem1 as rfem1,",
             "lfihk.beu_r_fem2 as rfem2,",
             "lfihk.beu_l_daq1 as ldaq1,",
             "lfihk.beu_r_daq1 as rdaq1",
             "from sci%d%d_weighted as sci" % (horn, rad),
             "join sci%d%d%d as sci_d0 using(pointingID)" % (horn,rad,0),
             "join sci%d%d%d as sci_d1 using(pointingID)" % (horn,rad,1),
             "join raw_gains%d%d as rg using(pointingID)"% (horn, rad),
             "join smoothed_gains%d%d as sg using(pointingID)" % (horn, rad),
             "join dipole_%d%d as dipole using(pointingID)" % (horn, rad),
             "join hfi_hk_temperatures as hfihk using(pointingID)",
             "join lfi_hk_temperatures as lfihk using(pointingID)",
             "join pointings as pointing_info using(pointingID)",
             "where vref > 0",
             "and rg.gain > 0",
             "and (sci.pointingID < 10911 or sci.pointingID > 10920)",
             "order by sci.pointingID" ]))
        all_list=[list(a) for a in c]
        for a in all_list:
            for i in range(len(a)):
                if a[i]=='':
                    a[i]=0
        ch_data = np.array([tuple(a) for a in all_list],dtype=ch_data_dtype)
        all_data[ch.tag] = ch_data
    return all_data    

def get_ucds2(db="/global/homes/p/peterm/ucds-dx10.db"):
    MS = 'MS'
    conn = sqlite3.connect(db)
    cursor=conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables=cursor.fetchall()
    dfall=pd.read_sql('select * from '+tables[0][0],conn)
    dfall=dfall.set_index('pointingID')
    for table in tables[1:]:
        df=pd.read_sql('select * from '+ table[0],conn)
        if 'pointingID' in df:
            df=df.set_index('pointingID')
            for col in df.columns:
                tcol=table[0]+'_'+col
                dfall[tcol]=df[col]
    return dfall
def get_digits(instring):
    """
    function to return digits in instring
    """
    digits=''
    for char in instring:
        if char in string.digits:
            digits=digits+char
    return digits
    

def get_ucds3(db="/global/homes/p/peterm/ucds-dx10.db"):
    """
    function to read ucds into a pandas hierarchical dataframe
    returns two frames: hk and science
    """
    conn = sqlite3.connect(db)
    cursor=conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables=cursor.fetchall()
    dfall=pd.read_sql('select * from '+tables[0][0],conn)
    n=len(dfall.pointingID)
    hk=resize(['hk'],n)
    dfallhk=dfall.set_index('pointingID')
    firstsci=True
    for table in tables[1:]:
        digits=get_digits(table[0])
        print digits
        df=pd.read_sql('select * from '+ table[0],conn)
        if 'pointingID' in df:
            n=len(df.pointingID)
            if len(digits)<3:
                hk=resize(['hk'],n)
                df=df.set_index([hk,'pointingID'])
                dfallhk=dfallhk.append(df)
            if len(digits)==3:
                if digits[2]=='0':
                    indx=resize([digits[:2]+'M'],n)
                if digits[2]=='1':
                    indx=resize([digits[:2]+'S'],n)
                df=df.set_index([indx,'pointingID'])
                dfallsci=dfallsci.append(df)
            if len(digits)==4:
                if digits[2]=='0':
                    if digits[3]=='0':
                        indx=resize([digits[:2]+'M0'],n)
                    if digits[3]=='1':
                        indx=resize([digits[:2]+'M1'],n)
                if digits[2]=='1':
                    if digits[3]=='0':
                        indx=resize([digits[:2]+'S0'],n)
                    if digits[3]=='1':
                        indx=resize([digits[:2]+'S1'],n)
                df=df.set_index([indx,'pointingID'])
                if not(firstsci):
                    dfallsci=dfallsci.append(df)                 
                if firstsci:
                    dfallsci=df
                    firstsci=False
    return dfallhk,dfallsci    
  
def spectrogram(indata,sampsperspec=1000,minfreq=None):
    """
    make simple spectrogram, uses nps
    """
    n=len(indata)
    nspec=np.int(n/sampsperspec)
    ns=nspec*sampsperspec
    indata=reshape(indata[:ns],(nspec,sampsperspec))
    sgram=[]
    for i in range(nspec):
        z=nps(indata[i,:],1.,minfreq=minfreq)
        zlen=len(z[0])
        sgram.append(z[1])
    sgram=np.array(sgram)
    sgram=np.reshape(sgram.flatten(),(nspec,zlen))
    return sgram
    
def spectrogramx(indata1,indata2,sampsperspec=1000,minfreq=None):
    """
    make simple cross spectrogram, uses xps
    """
    n=len(indata1)
    nspec=np.int(n/sampsperspec)
    ns=nspec*sampsperspec
    indata1=reshape(indata1[:ns],(nspec,sampsperspec))
    indata2=reshape(indata2[:ns],(nspec,sampsperspec))
    sgram=[]
    for i in range(nspec):
        z=xps(indata1[i,:],indata2[i,:],1.,minfreq=minfreq)
        zlen=len(z[0])
        sgram.append(z[1])
    sgram=np.array(sgram)
    sgram=np.reshape(sgram.flatten(),(nspec,zlen))
    return sgram

def syncpidplot(mapfig=1,linefig=2,od=None,ring=None,color='red',det='LFI28S',linewidth=3,ucds=None,coord=None):
    """
    function to plot a qr ring on a map (figure number mapfig=1) and vertical lines on a line plot. 
    if KW OD is given, convertes to min an max pid and plots those using qr, and a single line on the OD assumed x axis of the plot
    if KW ring is given, just plots QR for that ring and corresponding x axis (assumed PID) of line plot (linefig=2) Default detector
    is LFI28s, put specific one for more precision. set color and LW as well. If UCDS is not passed, it's read (slow)
    """
    if ring!= None:
        plt.figure(mapfig)
        qr.graph.plot_ring(ring+237,qr.det_opening_angles[det],color=color,linewidth=linewidth,label=str(ring),coord=coord)
        plt.draw
        plt.figure(linefig)
        plt.axvline(x=ring,color=color,linewidth=linewidth)
    if od!=None:
        rings=od2pid_range(od,ucds=ucds,key=det)
        plt.figure(mapfig)
        qr.graph.plot_ring(rings[0]+237,qr.det_opening_angles[det],color=color,linewidth=linewidth,label=str(od),coord=coord)
        qr.graph.plot_ring(rings[1]+237,qr.det_opening_angles[det],color=color,linewidth=linewidth,label='',coord=coord)
        plt.draw
        plt.figure(linefig)
        plt.axvline(x=od,color=color,linewidth=linewidth)        
        
def removeqr(mapfig=1,color='red'):
    plt.figure(mapfig)
    a=plt.gca()
    for li in a.lines:
        if li.get_color()==color:
            li.remove()
    plt.draw
    
def spectrogram_cds(dlist):
    """
    create spectrogram from output of get_cds_key_value. assumed list of OD minute sampled values (approx 1441/day)
    """
    sg=[]
    freqs=[]
    ods=[]
    od=91
    for d in dlist:
        z=nps(d,1./60.,minfreq=.00004)
        sg.append(z[1])
        freqs.append(z[0]*1000.)
        ods.append(np.zeros(len(z[1]))+od)
        od+=1
    sgout={}
    sgout['spectrogram']=sg
    sgout['frequency']=freqs
    sgout['od']=ods
    return sgout

    
def get_cds_key_value(key='fourk_temp',value='cernox_4k',cdsdir='/project/projectdirs/planck/data/mission/lfi_raw_sync/CDS/'):
    fls=glob(cdsdir+'*.sav')
    lenfls=[len(f) for f in fls]
    fls=np.array(fls)
    lenfls=np.array(lenfls)
    fls=fls[lenfls==72]
    fls.sort()
    d=[]
    for f in fls:
        print f
        cds=readsav(f)
        if key in cds.keys():
            d.append(cds[key][value][0].flatten())
    return d
        
def get_cds(cdsdir='/project/projectdirs/planck/data/mission/lfi_raw_sync/CDS/'):
    fls=glob(cdsdir+'*.sav')
    lenfls=[len(f) for f in fls]
    fls=np.array(fls)
    lenfls=np.array(lenfls)
    fls=fls[lenfls==72]
    fls.sort()
    cds=readsav(fls[0],python_dict=True)
    for f in fls:
        print f
        cd=readsav(f,python_dict=True)
        for key in cd.keys():
            cds[key]=np.hstack([cds[key],cd[key][0]])
    return cds
    

def nps(s, Fs,minfreq=None):
    """
    returns two vectors, frequencies and PSD
    PSD is in units^s/Hz
    """
    if minfreq != None:
        nfft=np.min([len(s),np.int(2.*Fs/minfreq)])
        nfft=2**(np.int(np.log2(nfft)))
    elif minfreq == None:
        nfft=len(s)
        nfft=2**(np.int(np.log2(nfft)))
    #Pxx, freqs = plt.psd(s, NFFT=nfft, Fs = Fs)
    Pxx, freqs = mlab.psd(s, NFFT=nfft, Fs = Fs)
    #we hate zero frequency
    freqs=freqs[1:]
    Pxx=Pxx[1:]
    return freqs, Pxx 
    
def xps(s1,s2, Fs,minfreq=None):
    """
    like nps, but for cross spectra: returns two vectors, frequencies and PSD
    PSD is in units^s/Hz
    """
    if minfreq != None:
        nfft=np.min([len(s1),np.int(2.*Fs/minfreq)])
        nfft=2**(np.int(np.log2(nfft)))
    elif minfreq == None:
        nfft=len(s1)
        nfft=2**(np.int(np.log2(nfft)))
    Pxx, freqs = mlab.csd(s1,s2, NFFT=nfft, Fs = Fs)
    #we hate zero frequency
    freqs=freqs[1:]
    Pxx=Pxx[1:]
    return freqs, Pxx 
    
def nps_old(s, Fs,minfreq=None):
    """
    returns two vectors, frequencies and PSD
    PSD is in units^s/Hz
    """
    if minfreq != None:
        nfft=np.min([len(s),np.int(2.*Fs/minfreq)])
        nfft=2**(np.int(np.log2(nfft)))
    elif minfreq == None:
        nfft=len(s)
        nfft=2**(np.int(np.log2(nfft)))
    Pxx, freqs = plt.psd(s, NFFT=nfft, Fs = Fs)
    #we hate zero frequency
    freqs=freqs[1:]
    Pxx=Pxx[1:]
    return freqs, Pxx 

def thicklegendlines(legendname,thick=3):
    lglines=legendname.get_lines()
    for line in lglines:
        line.set_linewidth(thick)
    plt.draw()
   
    
def rebin(a, range,ratio=True):
    """Downsizes a 2d array by averaging, new dimensions must be integral factors of original dimensions
Credit: Tim http://osdir.com/ml/python.numeric.general/2004-08/msg00076.html"""
    s=a.shape
    if len(s)==2:
        M, N = a.shape
        m=range[0]
        n=range[1]
        if ratio==True:
            ar = a.reshape((M/m,m,N/n,n))
            aout=np.mean(np.mean(ar, 2), 0) 
        if ratio==False:
            ar=a.reshape(m,M/m,n,N/n)
            aout=np.mean(np.mean(ar, 2), 0) 
    if len(s)==1:
        M=a.shape[0]
        m=range[0]
        if ratio==True:
            ar = a.reshape((M/m,m))
            aout=np.mean(ar, 1) 
        if ratio==False:
            ar=a.reshape(m,M/m)
            aout=np.mean(ar, 1) 
    return aout


def psd_function(freqs,wnlevel,fknee,alpha):
    """ 
    function to generate model PSD function (for fitting).
    Inputs freqs(array), wnlevel (value), Fknee (value), Alpha (value, for
    simple 1/f give alpha = 1)
    result will be power spectrum for all freqs in freqs, with
    wnlevel (input should be units/sqrt(Hz)) alpha will be for Power spectrum
    slope, not amplitude spectrum
    """
    nf=len(freqs)
    psdoutput=(wnlevel**2)*(1.+ (freqs/fknee)**(-alpha))
    return(psdoutput)

def getmapname(release='dx8',band='30',survey='1'):
    """
    function to point at lfi map files
    """
    dxlfipath='/project/projectdirs/planck/data/mission/DPC_maps/'+release+'/lfi/'
    fname=glob(dxlfipath+'LFI_'+band+'*_survey_'+survey+'.fits')
    return(fname)


def get_spin_phase(od):
    """
    function to open AHF file and read in spin phase. Use version v2 found in 
    /project/projectdirs/planck/data/mission/AHF_v2/<ODNUMBER>
    returns obt,pid,spin_phase
    """
    ahf_dir='/project/projectdirs/planck/data/mission/AHF_v2'
    ahf_file=glob(ahf_dir+'/*'+str(od)+'/*high*.fits')
    ahf=pyfits.open(ahf_file[0])
    obt=ahf[1].data['obt_spl']
    spin_phase=ahf[1].data['spin_phase']
    pid=ahf[1].data['pointing_id']
    ahf.close()
    #convert pid's to useable numbers- don't have official offset here yet
    pidfixed=[]
    for p in pid:
        pidfixed.append(np.int64(p[2:7]))
    pid=np.array(pidfixed)
    return(obt,pid,spin_phase)

def unwrap(angle):
    """
    function to unwrap an angular variable, used before interpolation
    assume degrees
    """
    ind=arange(len(angle))
    h=where(np.diff(angle) < -180.)
    hi=ind[h]
    nr=len(h)
    for i in hi:
        angle[i+1:]=angle[i+1:]+360.
    return(angle)
    
def sync_spinphase(obt_data,diff_data,obt_ahf,spinphase,pid=None):
    """
    function to interpolate spinphase to obt_data time base
    assume both are obt and spinphase is on 0-360
    we need to include diff_data in order to truncate it to 
    common times
    """
    #first step, truncate to common timeranges
    maxmin=np.max([np.min(obt_data),np.min(obt_ahf)])
    minmax=np.min([np.max(obt_data),np.max(obt_ahf)])
    diff_data=diff_data[(obt_data>maxmin) & (obt_data< minmax)]
    obt_data=obt_data[(obt_data>maxmin) & (obt_data< minmax)]
    spinphase=spinphase[(obt_ahf >maxmin) & (obt_ahf < minmax)]
    obt_ahf=obt_ahf[(obt_ahf >maxmin) & (obt_ahf < minmax)]
    if pid != None: 
        pid=pid[(obt_ahf >maxmin) & (obt_ahf < minmax)]
    
    #now unwrap angle
    spinphase=unwrap(spinphase)
    spinphaseinterp=np.interp(obt_data,obt_ahf,spinphase)
    spinphaseinterp=np.mod(spinphaseinterp,360.)
    pidinterp=None
    if pid != None:
        pidinterp=np.interp(obt_data,obt_ahf,pid)
        pidinterp=np.int32(pidinterp)
    return(obt_data,diff_data,spinphaseinterp,pidinterp)
    
    
def phasebin_pids(diff_data,spinphasei,pidi,azbins=1800):
    """
    function to bin by spin phase for all data in each PID
    assume already run sync_spinphase
    """
    azw=360./(azbins*2.)
    azs=(np.arange(azbins)+1.)*360./azbins
    pidlist=set(pidi)
    psmap=np.zeros((len(pidlist),azbins,3),dtype=np.float)
    for i,pid in enumerate(pidlist):
        psmap[i,:,0]=azs
        for j,az in enumerate(azs):
            d=diff_data[(pidi == pid) & (np.abs(spinphasei-az) < azw)]
            psmap[i,j,1]=np.mean(d)
            psmap[i,j,2]=np.std(d)/len(d)
    
    return(psmap)
    
def subtract_signal(diff_data,spinphasei,pidi,psmap):
    """
    function to subtract binned data stored in PSMap from diff_data
    """
    pidlist=set(pidi)
    diff_data_sr=np.zeros(len(diff_data),dtype=np.float)
    angwidth=np.abs(psmap[0,10,0]-psmap[0,9,0])/2.
    for i,pid in enumerate(pidlist):
        for j,ang in enumerate(psmap[i,:,0]):
            h=where((np.abs(spinphasei-ang) < angwidth) & (pidi == pid))
            if len(diff_data[h]) > 0: 
                diff_data_sr[h]=diff_data[h]-psmap[i,j,1]
    return(diff_data_sr)
         
def get_signal_removed_data(od,rca,diode,gmf=None):
    """
    script to combine the tools for reading data, ahf file and signal removing
    """
    obt,diff=get_lfi_diff(od,rca,diode,gmf=gmf)
    obt_ahf,pid,spinphase=get_spin_phase(od)
    obt_i,diff_i,spinphase_i,pid_i=sync_spinphase(obt,diff,obt_ahf,spinphase,pid)
    psmap=phasebin_pids(diff_i,spinphase_i,pid_i)
    diff_sigremoved=subtract_signal(diff_i,spinphase_i,pid_i,psmap)
    return(obt_i,diff_sigremoved)


    
def get_lfi_sky_ref(od,rca,diode,fits_dir='/project/projectdirs/planck/data/mission/lfi_ops_dx10'):
    """
    function to open fits file and read in sky and reference. data found in 
    /project/projectdirs/planck/data/mission/lfi_ops_dx10/
    rca  is number from 18 to 27
    diode is '00','01','10' or '11'  must be string!
    returns obt, sky, ref
    """
    #fits_dir='/project/projectdirs/planck/data/mission/lfi_raw_fits'
    fits_file=glob(fits_dir+'/od'+str(od)+'_rca'+str(rca)+str(diode)+'.fits')
    d=pyfits.open(fits_file[0])
    obt=d[1].data['time']
    sky=d[1].data['sky']
    ref=d[1].data['ref']
    d.close()
    return(obt,sky,ref)
    
def get_lfi_diff(od,rca,diode,gmf=None):
    """
    function to return obt and difference 
    data found in 
    /project/projectdirs/planck/data/mission/lfi_raw_fits/
    rca  is number from 18 to 27
    diode is '00','01','10' or '11'  must be string!
    returns obt, sky, ref
    """
    fits_dir='/project/projectdirs/planck/data/mission/lfi_raw_fits'
    fits_file=glob(fits_dir+'/od'+str(od)+'_rca'+str(rca)+str(diode)+'.fits')
    d=pyfits.open(fits_file[0])
    obt=d[1].data['time']
    sky=d[1].data['sky']
    ref=d[1].data['ref']
    d.close()
    if gmf==None:
        gmf=mean(sky)/mean(ref)
    diff=sky-gmf*ref
    return(obt,diff)
    
def get_lfi_dpc_timeline(od,rca,caltype='C',flag=False):
    """
    function to return obt and data from DPC timelines. 'C' means volts, use 'R' for reduced meaning Kcmb
    data found here: /project/projectdirs/planck/data/mission/lfi_ops_dx9_delta/XXXX (4 digit OD directory name)
    rca  is string of form 'LFI18M' if flag==true, skip flagged data
    returns obt, data
    """
    horn=int(rca[3:5])
    freq=30
    if horn<27:
        freq=44
    if horn<24:
        freq=70
        
    odstr=str(od)
    if len(odstr)==2:
        odstr='00'+odstr
    if len(odstr)==3:
        odstr='0'+odstr
    #fits_dir='/project/projectdirs/planck/data/mission/lfi_ops_dx9_delta/'+odstr
    fits_dir='/project/projectdirs/planck/data/mission/lfi_ops_dx10/'+odstr
    filelist=glob(fits_dir+'/L0'+str(freq)+'*'+caltype+'*.fits')
    filelist.sort()
    fits_file=filelist[-1]
    f=pyfits.open(fits_file)
    obt=f['OBT'].data['OBT']
    obtflag=f['OBT'].data['FLAG']
    d=f[rca].data[rca]
    f.close()
    if flag:
        obt=obt[obtflag==2]
        d=d[obtflag==2]
    return(obt,d)
    
    
def psd_fit_function(p,x):
    # Parameter values are passed in "p"
    # for PSD p=[wnlevel,fknee,alpha]
    # form is f(x)=(p[0]**2)*(1.+ (x/p[1])**(-p[2]))
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    # model = psd_function(x,p[0],p[1],p[2])
    model=(p[0]**2)*(1.+(x/p[1])**(-p[2]))
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return(model)
    
def psd_fit_function_resid(p,x,y,err):
    # Parameter values are passed in "p"
    # for PSD p=[wnlevel,fknee,alpha]
    # form is f(x)=(p[0]**2)*(1.+ (x/p[1])**(-p[2]))
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    # model = psd_function(x,p[0],p[1],p[2])
    model=(p[0]**2)*(1.+(x/p[1])**(-p[2]))
    status = 0
    return((y-psd_fit_function(p,x))/err)

def fit_fknee(psd,freqs):
    """
    function to call scipy optimize.leastsq to fit PSD function, assumed output of
    nps function above
    """
    if min(freqs) == 0:
        freqs=freqs[1:]
        psd=psd[1:]
    nfreq=len(freqs)
    topfreqs=np.where(freqs>.8*np.max(freqs))
    err=(np.zeros(nfreq,dtype=float)+np.std(np.sqrt(psd[topfreqs])))#*(1/freqs)
    p=np.array([np.sqrt(np.mean(psd[topfreqs])),.15,1.0])
    m=optimize.leastsq(psd_fit_function_resid,p,args=(freqs,psd,err),full_output=1)
    pfinal=m[0]
    print 'wnlevel',pfinal[0]
    print 'Fknee' ,pfinal[1]
    print 'alpha' ,pfinal[2]
    return(m)
    
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    (from cookbook)
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]    
    
    
def make_1overf(nsamples=1000000,wnlevel=1.0,fknee=.1,alpha=1.0,samprate=35.0):
    """
    function to produce TOI with specified 1/f spectrum. note wn level is in amplitude/sqrt(Hz), but
    alpha defined for PSD, not ASD.
    """
    freq=np.arange(2+nsamples/2,dtype=np.float32)*(samprate/nsamples)
    freq=freq[1:]
    psfilt=np.sqrt(psd_function(freq,wnlevel,fknee,alpha))
    toi=np.random.normal(0,np.sqrt(samprate/2.),nsamples)
    ztoi=np.fft.rfft(toi)
    fztoi=ztoi*psfilt
    toi_output=np.fft.irfft(fztoi)
    return(toi_output) 

if __name__ =='__main__':
    import matplotlib.pyplot as plt
    array_100x100= np.arange(100) * np.arange(100)[:,None]
    array_10x10 = rebin(array_100x100, (10, 10))
    plt.figure()
    plt.contourf(array_100x100)
    plt.contourf(array_10x10)
    plt.colorbar()
from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import cPickle
import numpy as np
from glob import glob
import pyfits


"""
models
LFI23M:
sigma=503.397 uK sqrt(s), alpha=-1.072, f_knee= 30.17mHz
LFI24S:
sigma=398.225 uK sqrt(s), alpha=-0.910, f_knee= 73.03mHz
"""
def logbin(freq,ps,nbins):
    """
    function to bin in logarithmic frequency intervals
    output freq and ps arrays will be nbins long
    """
    logfreqs=logspace(np.log2(freq[1]),np.log2(np.max(freq)),num=nbins,endpoint=True,base=2)
    
    outps=[ps[0]]
    outfreqs=[(freq[0]+freq[1])/2.]
    lastfreq=outfreqs[0]
    for logfreq in logfreqs:
        outps.append(np.mean(ps[np.logical_and(freq >= lastfreq,freq<logfreq)]))
        #outfreqs.append(np.mean(freq[np.logical_and(freq >= lastfreq,freq<=logfreq)]))
        outfreqs.append((logfreq+lastfreq)/2.)
        lastfreq=logfreq
    return np.array(outfreqs),np.array(outps)

width=8.8
#for width in [18., 12., 8.8]:
f_ny_23m=76.8/2.
f_ny_24s=45./2.
f_ny_27m=32.5/2. #all from timescales page on LFI wiki
sigmalist=[503.397,398.225]
alphalist=[-1.072,-0.910]
fklist=[30.17,73.03]
for ch,f_ny,sigma,alpha,fk in zip(['LFI23M','LFI24S'],[f_ny_23m,f_ny_24s],sigmalist,alphalist,fklist):
    fl=glob('noise*'+ch+'.fits')
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    
    for f in fl:
        ff=pyfits.open(f)
        logps=logbin(1000.*ff[1].data['FREQ'],ff[1].data['PS'],20)
        plt.plot(logps[0],1000.*np.sqrt(logps[1]/f_ny))
        ff.close()
    x=logspace(np.log2(logps[0][0]),np.log2(logps[0][-1]),num=200,base=2)
    model=.001*np.sqrt(2*(sigma**2)*(1+(x/fk)**alpha))
    plt.plot(x,model,lw=2,color='black',label='Noise model')
    plt.xscale('log')
    plt.yscale('log')
        
    plt.legend(frameon=False,loc=0)
    plt.grid()
    
    plt.ylim([.3,10.])
    #plt.xlim([.0005,40])
    plt.show()
    # labels
    plt.xlabel(r"Frequency, mHz"); plt.ylabel(r"Amplitude Spectral Density $\left[ mK_{cmb} /\sqrt{Hz} \right]$")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # reduce white space around figure
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")
    plt.show()
    plt.savefig("asd_sample"+ch+"_%dmm.eps" % int(width*10), bbox_inches='tight',pad_inches=0)
    plt.savefig("asd_sample"+ch+"_%dmm.pdf" % int(width*10), bbox_inches='tight',pad_inches=0)
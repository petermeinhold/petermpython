"""set of utility functions for acquisition software"""

import numpy as np
import pyfits as pyfits
import ephem 
#import peakanalysis as pk
import matplotlib.pyplot as plt
#from planck_util import nps
from glob import glob
from scipy.signal import butter,filtfilt,iirdesign
#from demod import datparsing
from matplotlib import pyplot as plt

rtd=180/np.pi

def bcd_to_int(int_1,int_2):
	#function to pull out bytes from int_1 and int_2 and reconstruct BCD decimal number
	b_1=[]
	b_2=[]
	for i1 in int_1:
		b_1.append(bin(i1))
	for i2 in int_2:
		b_2.append(bin(i2))
	outinteger=[]
	bb1=0
	bb2=0
	bb3=0
	bb4=0
	for b1,b2 in zip(b_1,b_2):
		if len(b1)==10:
			bb1=int(b1[6:10],2)
			bb2=int(b1[2:6],2)
		if len(b1)==6:
			bb1=int(b1[2:6],2)
			bb2=0
		if len(b2)==10:
			bb3=int(b2[6:10],2)
			bb4=int(b2[2:6],2)
		if len(b2)==6:
			bb3=int(b2[2:6],2)
			bb4=0
		print bb1,bb2,bb3,bb4
		outinteger.append( bb1+ bb2*10+100*bb3+100*bb4)

	return outinteger
	
		
def phasebin(nbins, az, signal):
    ring_edges=np.where(np.diff(az) < -np.pi)
    nrings=len(ring_edges[0])
    phasebin_edges = np.linspace(-np.pi,np.pi, nbins+1)
    #SPLIT THE PHASE FOR EACH RING
    pseudomap = np.zeros([nbins,nrings-1],dtype=np.float32 ) 
    #plt.figure()
    #plt.hold(False)
    for ring in range(nrings-1):
        az_ring=az[ring_edges[0][ring]:ring_edges[0][ring+1]]
        signal_ring=signal[ring_edges[0][ring]:ring_edges[0][ring+1]]
        pseudomap[:,ring], edges = np.histogram(az_ring, bins=phasebin_edges, weights=signal_ring)
        #plt.plot(edges[:-1],pseudomap[:,ring])
        #plt.show()
        hits, edges = np.histogram(az_ring, bins=phasebin_edges)
        pseudomap[hits>0,ring] /= hits[hits>0]
    return pseudomap

def thicklegendlines(legendname,thick=3):
    lglines=legendname.get_lines()
    for line in lglines:
        line.set_linewidth(thick)
    plt.draw()
    
def lowpass(d,sample_rate,cutoff):
    '''
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html#scipy.signal.filtfilt
    just stuff in an example from scipy to get this functional
    '''
    frac_cutoff=cutoff/(sample_rate/2.)
    print frac_cutoff
    b,a=butter(3,frac_cutoff)
    #b,a=iirdesign(frac_cutoff-.001,frac_cutoff+.1,.9,.1)
    filtered_d = filtfilt(b,a,d)
    return(filtered_d)

def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))

def rebin_factor( a, newshape ):
        '''Rebin an array to a new shape.
        newshape must be a factor of a.shape.
        '''
        assert len(a.shape) == len(newshape)
        assert not np.sometrue(np.mod( a.shape, newshape ))

        slices = [ slice(None,None, old/new) for old,new in zip(a.shape,newshape) ]
        return a[slices]

if __name__ =='__main__':
    import matplotlib.pyplot as plt
    array_100x100= np.arange(100) * np.arange(100)[:,None]
    array_10x10 = rebin(array_100x100, (10, 10))
    plt.figure()
    plt.contourf(array_100x100)
    plt.contourf(array_10x10)
    plt.colorbar()
    
def remove_bad_raw_data(d):
    """
    function to find noise triggered data 
    returns cut down memmap whie still needs to be0
    tested for incomplete revolutions etc.
    """
    encoder=d['enc']/16
    g=np.where(np.diff(encoder) != 0)
    return (d[g])

def get_raw_data(utmin,utmax,freq='10'):
    '''
    function to turn the UT range into filenames and read in the raw data
    '''
    ltmin=utmin-7.0 #convert to local time first
    ltmax=utmax-7.0
    data_dir='/cofe/flight_data/'+str(freq)+'GHz/'
    if (ltmin<24) & (ltmax <24):
        data_dir=data_dir+'20110917/'     
    elif (ltmin>24) & (ltmax >24):
        data_dir=data_dir+'20110917/'
        ltmin=ltmin-24.
        ltmax=ltmax-24.
    
    #list the available files
    fl=glob(data_dir+'*.dat')
    ltfl=[]
    for file in fl:
        ltfl.append(file[-12:-4])
    ltflhours=np.zeros(len(ltfl))
    for i,lt in enumerate(ltfl):
        ltflhours[i]=np.float(lt[0:2])+np.float(lt[2:4])/60.+np.float(lt[4:6])/3600.
    fl=np.array(fl)
    ltflhours=np.array(ltflhours,dtype=float)
    files2read=fl[(ltflhours>ltmin) & (ltflhours<ltmax)]
    len(files2read)
    d=datparsing.read_raw(files2read)
    return(d)
    
    
    
def get_cofe_target(ut,lat,lon,target):
#function to use ephem to find az and el of specified target for COFE
#parameters: UT, LAt Lon Target
    cofe=ephem.Observer()
    cofe.elevation=300000
    year=2011
    month=9
    day=17
    az=[]
    el=[]
    for u,la,lo in zip(ut,lat,lon):
        if (u >24):
            u=u-24
            day=18
        hour=int(np.fix(u))
        minute=(u-hour)*60
        iminute=int(np.fix(minute))
        second=int(np.fix((minute-iminute)*60))
        datestring=str(year)+'/'+str(month)+'/'+str(day)+' '+str(hour)+':'+str(iminute)+':'+str(second)
        datestring
        cofe.date=datestring
        cofe.lon=str(rtd*lo)
        cofe.lat=str(rtd*la)
        pos=ephem.__getattribute__(target)(cofe)
        az.append(pos.az)
        el.append(pos.alt)
        
    return np.array(az),np.array(el)
    
def get_cofe_crossing(ut,toi,gaz,lat,lon,centerut,target,plot=False):
    #function to find maximum peak signal crossing
    #in toi, return azoffset, elevation at crossing
    t=np.where(abs(ut-centerut) < .3/60.)
    h,=np.where(toi[t] == np.min(toi[t]))
    h=np.int(np.median(h))
    targetpos=get_cofe_target(ut[t],lat[t],lon[t],target)
    gfit=pk.fit_a_peak(toi[t],h,invert=1)
    print gfit[0].params
    x=np.arange(1,len(ut[t])+1)
    gauss=pk.gaussian(gfit[0].params,x)
    if plot==True:
        plt.plot(x,toi[t])
        plt.hold(True)
        plt.plot(x,gauss)
    center=np.round(gfit[0].params[2])
    targetaz=rtd*targetpos[0][center]
    targetel=rtd*targetpos[1][center]
    azoffset=targetaz-gaz[center+t[0][0]]*rtd
    fitut=ut[center+t[0][0]]
    return azoffset,targetaz,targetel,fitut
    
def cctout(cc, gpstime):
    """cc and gpstime must be already synched"""
    utc = np.mod(((gpstime+15.)/3600.), 24)
    ut = utc[0]+(((cc - cc[0])*2e-9)/3600.)
    return ut
    
    
def oplot(*params):
    """function to emulate one good feature of IDL oplot with nice features of Python plot"""
    plt.hold(True)
    plot(*params)
    plt.hold(False)
    return
    
    
def linfit(x,y):
    """embed python lin algebra fitting to look like idl tool"""
    if len(y) != len(x):
        print 'inputs need to be same length arrays'
        
    a=np.vstack([x,np.ones(len(x))]).T
    m,b=np.linalg.lstsq(a,y)[0]
    return np.array([b,m])

def find_command_uts(command):
#script to read in cmdecho file, find lines with 
#requested command, and return adjusted UT (adjusted
#by offset to match v 0.5 level 1 files, calibration signal
# in ch a/d 1 of 15 GHz. near 16:42 UT.
    import datetime
    f=open('C:\\cofe\\flight_data\\CIP_commanding\\623nmisc\\cmdecho.dat')
    
    timelist=[]
   
    for line in f:
        linelist=line.split("\t")
        linedate=datetime.datetime.strptime(linelist[0]+' '+linelist[1],'%m/%d/%Y %H:%M:%S')
        addr=linelist[2]
        cmd=(linelist[3].strip('\n'))
        if ((cmd == command) and (addr == '0006')):
            print 'found one'
            #utadjusted=24*(linedate.day-17.)+linedate.hour+linedate.minute/60.+linedate.second/3600.
            #timelist.append(utadjusted)
            timelist.append((linedate.day,linedate.hour,linedate.minute,linedate.second))
    return timelist
    f.close()
    
def raise_bit(flag, bit=0):
    '''Raise bit of the flag array'''
    return np.bitwise_or(flag, 2**bit)

def check_bit(flag, bit=0):
    '''Check if bit of the flag array is raised'''
    return np.bitwise_and(flag, int(2**bit)) > 0
    
def grab_x_from_plot(fig):

    print 'right button press selects xvalue to store. middle click to  end function'
    global startlist,stoplist,start
    global cid
    global ptnum
    ptnum=0
    start=True
    startlist=[]
    stoplist=[]
    
    def onclick(event):
        global ptnum,startlist,stoplist,start
        global cid
        print ptnum
        print ptnum%2
        if event.button == 3:
            if ptnum%2 == 0:
                startlist.append(event.xdata)
                print 'chose start point',event.xdata
                print 'select stop point'
            elif ptnum%2 == 1:
                stoplist.append(event.xdata)
                print 'chose stop point',event.xdata
                print 'select next startpoint'
            start=False
            ptnum+=1
        elif event.button ==2:
            print 'should quit now, was here: ',event.xdata
            print cid
            fig.canvas.mpl_disconnect(cid)
    cid=fig.canvas.mpl_connect('button_press_event',onclick)
    outlist=np.concatenate((np.array(startlist),np.array(stoplist)),axis=2)
    return(outlist)

            
    
    
    
    
    
    
        
        
    
        
        

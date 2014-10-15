
import numpy as np
import pyfits as pyfits
import ephem 
import matplotlib.pyplot as plt
from glob import glob
from matplotlib import pyplot as plt
from matplotlib import mlab
import time
#import quaternion as quat
import numexpr
import numpy

# Define a local second in Julian time - useful for duration calculation
second_jul = 0.1 / 8640

# Define pi - useful for evaluation expressions
pi = numpy.pi
rtd=180/np.pi
dtr=np.pi/180

def eqmap2azelmap(map,longi=-104.245,lat=34.4717,ut=12.0,year=2014,month=9,day=15):
    """
    function to rotate a galactic coords healpix map to azel coordinates given
    longitude,latitude,altitude and time. defaults are 30 km above fort sumner
    Idea is to use pyephem to find the galactic coords of local zenith, then use Andrea's quaternion
    routines to find the rotation matrix, then use this to rotate all of the heapix pixels to azel.
    """
    nside=np.int(np.sqrt(len(map)/12))
    obs=ephem.Observer()
    obs.lon=str(longi)
    obs.lat=str(lat)
    obs.date=ephem.date((year,month,day,ut))
    radeczenith=ephem.Equatorial(*obs.radec_of(0,np.pi/2.))   #'unpack' the tuple produced by obs.radec_of.. result here is pyephem object class
    theta=np.pi/2.-radeczenith.dec
    phi=radeczenith.ra
    veczenith=hp.ang2vec(theta,phi) #should be healpix convention vector of gal coords of local zenith
    vecaxis=hp.ang2vec(0,np.pi)                                       #ref line in gal coords (ie rotation matrix should take veczenith here)    
    q=quat.from_vectors(veczenith,vecaxis)
    rmat=quat.to_rotmat(q)
    pixlist=range(len(map))
    veclist=hp.pix2vec(nside,pixlist)
    rveclist=np.dot(rmat,veclist)
    rpixlist=hp.vec2pix(nside,rveclist[0,:],rveclist[1,:],rveclist[2,:])
    outmap=zeros(len(map))
    outmap[rpixlist]=map[pixlist]
    return outmap
    
    
def eqmap2azelmap_brute(map,longi=-104.245,lat=34.4717,ut=12.0,year=2014,month=9,day=15):
    """
    function to rotate celestial coord map to az el with brute force using tools below from Victor Roytman
    """
    julian_date=jdcnv(year,month,day,ut)
    obs=ephem.Observer()
    obs.lon=str(longi)
    obs.lat=str(lat)
    obs.date=ephem.date((year,month,day,ut))
    nside=np.int(np.sqrt(len(map)/12))
    outmap=np.zeros(len(map))
    pixlist=range(len(map))
    htheta,hphi=hp.pix2ang(nside,pixlist)
    elevation=np.pi/2. -htheta
    azimuth=hphi
    ctheta=[]
    cphi=[]
    for az,el in zip(azimuth,elevation):
        ra,dec=obs.radec_of(az,el)
        ctheta.append(np.pi/2. -dec)
        cphi.append(ra)
    #ra,dec=azel2radec(julian_date, azimuth, elevation, lat*dtr, longi)
    #ctheta=np.pi/2.-dec
    #cphi=ra
    ctheta=np.array(ctheta)
    cphi=np.array(cphi)
    rpixlist=hp.ang2pix(nside,ctheta,cphi)
    print len(rpixlist)
    outmap[pixlist]=map[rpixlist]
    return outmap
    
def galmap2eqmap(map):
    """
    function to rotate galactic coord map (healpix) to equatorial
    """
    nside=np.int(np.sqrt(len(map)/12))
    grot=hp.Rotator(coord='GC')
    pixlist=range(len(map))
    veclist=hp.pix2vec(nside,pixlist)
    rveclist=grot(veclist)
    rpixlist=hp.vec2pix(nside,rveclist[0,:],rveclist[1,:],rveclist[2,:])
    outmap=zeros(len(map))
    outmap[rpixlist]=map[pixlist]
    return outmap
    

def jdcnv(year, month, day, hour):
    """Convert Gregorian time (UTC) to Julian Date
    Input:
    year = year (scalar int)
    month = month 1-12 (scalar int)
    day = day 1-31 (scalar int)
    hour = fractional hour (scalar double)
    Output:
    julian = Julian Date (scalar double)
    Original IDL source available at:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/jdcnv.pro
    Revision history of IDL source:
    Converted to IDL from Don Yeomans Comet Ephemeris Generator,
    B. Pfarr, STX, 6/15/88
    Converted to IDL V5.0 W. Landsman September 1997
    Added checks on valid month, day ranges W. Landsman July 2008
    """
    year = long(year)
    month = long(month)
    day = long(day)
    # account for leap years
    leap = long((month - 14) / 12.0)
    
    julian = day - 32075L + long((1461L) * (year + 4800L + leap) / 4.0) \
             + long((367L) * (month - 2 - leap * 12) / 12.0) \
             - long(3 * (long((year + 4900L + leap) / 100.0)) / 4.0) \
             + (hour / 24.0) - 0.5
    
    return julian

def ct2lst(julian_date):
    """Convert Civil Time (as Julian date) to Greenwich Sidereal Time
    Input:
    julian_date = Julian date (scalar or vector double)
    Output:
    gst = Greenwich Sidereal Time (scalar or vector double)
    The constants used in ct2lst come from Astronomical Algorithms by Jean
    Meeus, p. 84 (Eq. 11-4).
    Original IDL source available at:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro
    Revision history of IDL source:
    Adapted from the FORTRAN program GETSD by Michael R. Greason, STX,
    27 October 1988.
    Use IAU 1984 constants Wayne Landsman, HSTX, April 1995, results
    differ by about 0.1 seconds
    Longitudes measured *east* of Greenwich W. Landsman December 1998
    Time zone now measure positive East of Greenwich W. Landsman July 2008
    Remove debugging print statement W. Landsman April 2009
    """
    c1 = 280.46061837
    c2 = 360.98564736629
    c3 = 0.000387933
    c4 = 38710000.0
    t0 = numexpr.evaluate('julian_date - 2451545.0')
    t = numexpr.evaluate('t0 / 36525')
    
    # Compute GST in seconds
    theta = numexpr.evaluate('c1 + (c2 * t0) + (t**2)'\
                             ' * (c3 - (t / c4))')
    
    # Compute LST in hours
    gst = numexpr.evaluate('theta / 15.0')
    
    # Deal with LST out of bounds
    negative = numexpr.evaluate('gst < 0.0')
    negative_gst = gst[negative]
    negative_gst = numexpr.evaluate('24.0 + negative_gst % 24')
    gst[negative] = negative_gst
    gst = numexpr.evaluate('gst % 24.0')
    
    return gst

    # IMPORTANT: Latitude is in RADIANS
    # Longitude is in DEGREES
    # in order to reduce function calls
def azel2radec(julian_date, azimuth, elevation, latitude, longitude_deg):
    """Convert from horizon coordinates (azimuth-elevation) to celestial
    coordinates (right ascension-declination)
    Input:
    julian_date = Julian date (scalar or vector double)
    azimuth = azimuth in radians (scalar or vector double)
    elevation = elevation in radians (scalar or vector double)
    latitude = latitude in radians north of the equator (scalar or
    vector double)
    longitude_deg = longitude in degrees east of the prime meridian
    (scalar or vector double)
    Output:
    ra = right ascension in radians (scalar or vector double)
    dec = declination in radians (scalar or vector double)
    Original IDL source available at:
    http://cosmology.berkeley.edu/group/cmbanalysis/forecast/idl/azel2radec.pro
    Revision history of IDL source:
    Created, Amedeo Balbi, August 1998 (based partly on material
    by Pedro Gil Ferreira)
    Modified, Amedeo Balbi, October 1998, to accept vectors as input
    """
    
    # Rescale azimuth
    azimuth = numexpr.evaluate('azimuth - pi')
    
    # Get the Greenwich Sidereal Time of the date
    gst = ct2lst(julian_date)
    
    # Get the Local Sidereal Time
    lst = numexpr.evaluate('gst + longitude_deg / 15.0')
    
    # Calculate the declination
    dec = numexpr.evaluate('arcsin( '\
                           'sin(elevation) * sin(latitude) - '\
                           'cos(elevation) * cos(azimuth) * cos(latitude)'\
                           ' )')
    
    # Calculate the right ascension from the hour angle
    ha = numexpr.evaluate('arctan2( '\
                          '-(cos(elevation)*sin(azimuth)) , '\
                          '(sin(elevation)*cos(latitude)-'\
                          'cos(elevation)*cos(azimuth)*sin(latitude))'\
                          ' )')
    ha = numexpr.evaluate('24 * ha / 2.0 / pi')
    ra = numexpr.evaluate('lst - ha')
    ra=lst*2*pi/24. - ha
    #Deal with RA out of bounds
    negative = numexpr.evaluate('ra<0.0')
    negative_ra = ra[negative]
    negative_ra = numexpr.evaluate('24.0 + negative_ra % 24')
    ra[negative] = negative_ra
    ra = numexpr.evaluate('ra % 24.0')
    
    # Convert from hours to radians
    ra = numexpr.evaluate('ra * pi / 12')
    
    return ra, dec
    
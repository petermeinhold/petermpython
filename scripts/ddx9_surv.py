from petermpython import planck_util as pu
import numpy as np
import cPickle
import pycfitsio
from collections import defaultdict
import healpy as hp
from planck import Planck
from glob import glob

base='/global/project/projectdirs/planck/user/zonca/issues/dx9_toast/maps/'   
fl70a=glob(base+'*map*70*k_1*.fits')
fl70b=glob(base+'*map*70*k_2*.fits')
fl70a.sort
fl70b.sort

m70s1at=hp.ma(hp.read_map(fl70a[0],0))
m70s1aq=hp.ma(hp.read_map(fl70a[0],1))
m70s1au=hp.ma(hp.read_map(fl70a[0],2))

m70s2at=hp.ma(hp.read_map(fl70a[1],0))
m70s2aq=hp.ma(hp.read_map(fl70a[1],1))
m70s2au=hp.ma(hp.read_map(fl70a[1],2))

m70s3at=hp.ma(hp.read_map(fl70a[2],0))
m70s3aq=hp.ma(hp.read_map(fl70a[2],1))
m70s3au=hp.ma(hp.read_map(fl70a[2],2))

m70s4at=hp.ma(hp.read_map(fl70a[3],0))
m70s4aq=hp.ma(hp.read_map(fl70a[3],1))
m70s4au=hp.ma(hp.read_map(fl70a[3],2))

tmaska=m70s1aq.mask|m70s2aq.mask|m70s3aq.mask|m70s4aq.mask


m70s1bt=hp.ma(hp.read_map(fl70b[0],0))
m70s1bq=hp.ma(hp.read_map(fl70b[0],1))
m70s1bu=hp.ma(hp.read_map(fl70b[0],2))

m70s2bt=hp.ma(hp.read_map(fl70b[1],0))
m70s2bq=hp.ma(hp.read_map(fl70b[1],1))
m70s2bu=hp.ma(hp.read_map(fl70b[1],2))

m70s3bt=hp.ma(hp.read_map(fl70b[2],0))
m70s3bq=hp.ma(hp.read_map(fl70b[2],1))
m70s3bu=hp.ma(hp.read_map(fl70b[2],2))

m70s4bt=hp.ma(hp.read_map(fl70b[3],0))
m70s4bq=hp.ma(hp.read_map(fl70b[3],1))
m70s4bu=hp.ma(hp.read_map(fl70b[3],2))

tmaskb=m70s1bq.mask|m70s2bq.mask|m70s3bq.mask|m70s4bq.mask


nside=1024
freq=70
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/' + 'mask_ps_%dGHz_*.fits' % freq)[0]), nside)))

tmask=tmaska|tmaskb|psmask|wmask

y1at=(m70s1at+m70s2at)/2.0
y1bt=(m70s1bt+m70s2bt)/2.0
y2at=(m70s3at+m70s4at)/2.0
y2bt=(m70s3bt+m70s4bt)/2.0
y1at.mask=tmask
y1bt.mask=tmask
y2at.mask=tmask
y2bt.mask=tmask

y1aq=(m70s1aq+m70s2aq)/2.0
y1bq=(m70s1bq+m70s2bq)/2.0
y2aq=(m70s3aq+m70s4aq)/2.0
y2bq=(m70s3bq+m70s4bq)/2.0
y1aq.mask=tmask
y1bq.mask=tmask
y2aq.mask=tmask
y2bq.mask=tmask

y1au=(m70s1au+m70s2au)/2.0
y1bu=(m70s1bu+m70s2bu)/2.0
y2au=(m70s3au+m70s4au)/2.0
y2bu=(m70s3bu+m70s4bu)/2.0
yat=(y1at+y2at)/2.0
yaq=(y1aq+y2aq)/2.0
yau=(y1au+y2au)/2.0

ybt=(y1bt+y2bt)/2.0
ybq=(y1bq+y2bq)/2.0
ybu=(y1bu+y2bu)/2.0

ya=[yat,yaq,yau]
yb=[ybt,ybq,ybu]


y1abt=(y1at+y1bt)/2.0
y2abt=(y2at+y2bt)/2.0
y1abq=(y1aq+y1bq)/2.0
y2abq=(y2aq+y2bq)/2.0
y1abu=(y1au+y1bu)/2.0
y2abu=(y2au+y2bu)/2.0
y1ab=[y1abt,y1abq,y1abu]
y2ab=[y2abt,y2abq,y2abu]
cly12tt,cly12ee,cly12bb,a,b,c=hp.anafast(y1ab,y2ab)


y1au.mask=tmask
y1bu.mask=tmask
y2au.mask=tmask
y2bu.mask=tmask

y1a=[y1at,y1aq,y1au]
y1b=[y1bt,y1bq,y1bu]
y2a=[y2at,y2aq,y2au]
y2b=[y2bt,y2bq,y2bu]

cly1abtt,cly1abee,cly1abbb,a,b,c=hp.anafast(y1a,y1b)
cly2abtt,cly2abee,cly2abbb,a,b,c=hp.anafast(y2a,y2b)
cly12att,cly12aee,cly12abb,a,b,c=hp.anafast(y1a,y2a)
cly12btt,cly12bee,cly12bbb,a,b,c=hp.anafast(y1b,y2b)

l=arange(len(cly1abtt))
frac=float(len(y1au[tmask==False]))/float(len(tmask))

figure()
plot l,l*(l+1)*cly1abee/frac,label='Y1 half ring AxB EE'
plot l,l*(l+1)*cly2abee/frac,label='Y2 half ring AxB EE'
plot l,l*(l+1)*cly12aee/frac,label='Y1xY2 ringhalf1 EE'
plot l,l*(l+1)*cly12bee/frac,label='Y1 X Y2 ringhalf2 EE'

leg=legend()
pu.thicklegendlines(leg)
title('DX9 70 GHz Xspectra ),xlabel('Multipole'),ylabel('l*(l+1)Cl/2pi, microK^2')      

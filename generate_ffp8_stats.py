import sys
sys.path.append('/global/homes/p/peterm/petermpython/paperplots/python/scripts')
from setup_matplotlib import *
import healpy as hp
import planck_util_prm as pu
#get mask, assume lfi maps at nside 1024
nside=1024
tmask=pu.get_lfi_dx11_mask(nside)

#Create stats
topdir='/global/project/projectdirs/planck/data/ffp8/mc_noise/'
freqs=['070','044','030']
quads=['_18_23','_19_22','_20_21']
s1list=['full','full','full','full','full','y13']
s2list=['full','y1','y2','y3','y4','y24']
hr1list=['_hr1','','','','','']
hr2list=['_hr2','','','','','']
s1list=['full','full','full','full','full','full','full','full']
s2list=['s1','s2','s3','s4','s5','s6','s7','s8']
hr1list=['','','','','','','','']
hr2list=['','','','','','','','']
q1list=['','','','','','','','']
q2list=['','','','','','','','']

s1list=['full','full','full','full','full','full','y1','y1','y1','y4','y4','y4']
s2list=['y4','y4','y4','full','full','full','y1','y1','y1','y4','y4','y4']
hr1list=['','','','','','','','','','','','']
hr2list=['','','','','','','','','','','','']
q1list=['_18_23','_19_22','_20_21','_18_23','_18_23','_19_22','_18_23','_18_23','_19_22','_18_23','_18_23','_19_22']
q2list=['_18_23','_19_22','_20_21','_19_22','_20_21','_20_21','_19_22','_20_21','_20_21','_19_22','_20_21','_20_21']
freq='070'
mcdirlist=sort(glob(topdir+freq+'/'+freq+'_0*'))
homedir='/global/homes/p/peterm/'
for q1,q2,s1,s2,hr1,hr2 in zip(q1list,q2list,s1list,s2list,hr1list,hr2list):
    cls=[]
    for mcdir in [mcdirlist[0]]:
        #extract the 'hundreds' column to build list of mc numbers in this directory
        h=mcdir[-1]
        nstart=100*int(h)
        mcnumlist=[str(i).zfill(5) for i in range(nstart,nstart+100)]
        for mcnum in mcnumlist:
            f1=mcdir+'/ffp8_noise_%s%s_%s_map_mc_%s%s.fits' %(freq,q1,s1,mcnum,hr1)
            f2=mcdir+'/ffp8_noise_%s%s_%s_map_mc_%s%s.fits' %(freq,q2,s2,mcnum,hr2)
            cls.append(read_and_diff_files(f1,f2,nside=256,tmask=tmask))
    pklfilename=homedir+'ffp8_noise_cls/maps/'+'ffp8_noise_null_cls_'+freq+q1+s1+hr1+q2+s2+hr2+'.pkl'
    pklfile=open(pklfilename,'wb')
    cPickle.dump(cls,pklfile)
    pklfile.close()
    ffpcl_means=fitffpcls(cls)
    outfilename=pklfilename.replace('ffp8_noise_null_cls_','cl_fit_ffp8_noise_null_cls_')
    pklfile=open(outfilename,'wb')
    cPickle.dump(ffpcl_means,pklfile)
    pklfile.close()

from glob import glob
import healpy as hp
base='/global/homes/z/zonca/m/maps/out'
fl=glob(base+'/map*e1*.fits')
m1=hp.ma(hp.read_map(fl[7],[0,1,2]))
m2=hp.ma(hp.read_map(fl[9],[0,1,2]))
nside=1024
freq=70
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
tmask=wmask.astype(np.bool)|psmask|m1[1].mask|m2[1].mask
for m in m1:
    m.mask=tmask
for m in m2:
    m.mask=tmask
cl_e1=hp.anafast(m1,m2)
l=np.arange(len(cl_e1[0]))
frac=(1-float(m1[0].mask.sum())/len(m1[0].mask))

fla=glob(base+'/map*e3*.fits')
m1a=hp.ma(hp.read_map(fla[2],[0,1,2]))
m2a=hp.ma(hp.read_map(fla[0],[0,1,2]))
tmaska=wmask.astype(np.bool)|psmask|m1a[1].mask|m2a[1].mask
for m in m2a:
    m.mask=tmaska
for m in m1a:
    m.mask=tmaska
cl_e3=hp.anafast(m1a,m2a)
fldx9=glob(base+'/map_ddx9*18_23*full*subchunk*.fits')
m1dx9=hp.ma(hp.read_map(fldx9[1],[0,1,2]))
m2dx9=hp.ma(hp.read_map(fldx9[0],[0,1,2]))
tmaskdx9=wmask.astype(np.bool)|psmask|m1dx9[1].mask|m2dx9[1].mask
for m in m1dx9:
    m.mask=tmaskdx9
for m in m2dx9:
    m.mask=tmaskdx9
cl_dx9=hp.anafast(m1dx9,m2dx9)
figure()
plot(l,1e12*cl_e1[4]*l*(l+1)/(2*np.pi*frac),label='EB Mode e1')
plot(l,1e12*cl_e3[4]*l*(l+1)/(2*np.pi*frac),label='EB Mode e3')
leg=legend()
plot(l,1e12*cl_dx9[4]*l*(l+1)/(2*np.pi*frac),label='EB DDX9')
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EB halfRing cross-spectra for ELina cal mode 1, mode 2 and DDX9, Full')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-20,20])
plt.savefig('cl_elina_cal_EB.png')

plt.figure()
plt.plot(l,1e12*cl_e1[1]*l*(l+1)/(2*np.pi*frac),label='EE Mode e1')
plt.plot(l,1e12*cl_e3[1]*l*(l+1)/(2*np.pi*frac),label='EE Mode e3')
leg=plt.legend()
plt.plot(l,1e12*cl_dx9[1]*l*(l+1)/(2*np.pi*frac),label='EE DDX9')
plt.xlabel('Multipole'),plt.ylabel('l(l+1)Cl, microK^2'),plt.title('EE halfRing cross-spectra for ELina cal mode 1, mode 2 and DDX9, Full')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([0,20])
plt.savefig('cl_elina_cal_EE.png')

figure()
plot(l,1e12*cl_e1[2]*l*(l+1)/(2*np.pi*frac),label='BB Mode e1')
plot(l,1e12*cl_e3[2]*l*(l+1)/(2*np.pi*frac),label='BB Mode e3')
leg=legend()
plot(l,1e12*cl_dx9[2]*l*(l+1)/(2*np.pi*frac),label='BB DDX9')
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('BB halfRing cross-spectra for ELina cal mode 1, mode 2 and DDX9, Full')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([0,20])
plt.savefig('cl_elina_cal_BB.png')
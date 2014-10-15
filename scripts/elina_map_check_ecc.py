
from glob import glob
import healpy as hp
base='/global/project/projectdirs/planck/user/zonca/issues/dx9_delta_toast/maps/out/'
surveys=['survey1','survey2','survey3','survey4','survey5']
nside=1024
freq=70
psmask = np.logical_not(np.floor(hp.ud_grade(hp.read_map(glob('/project/projectdirs/planck/data/mission/DPC_maps/dx8/lfi/DX8_MASKs/mask_ps_%dGHz_*.fits' % freq)[0]), nside)))
wmask =(hp.ud_grade(hp.read_map('/project/projectdirs/planck/user/zonca/issues/singler/wmap_polarization_analysis_mask_r9_7yr_v4.fits'),nside) < 1).astype(np.bool)
#'ecc' or 'ddx9_18'
figure(4)
figure(5)
figure(6)
for survey in surveys:
    fl1=glob(base+'map_ecc*'+survey+'_subchunk_1*.fits')
    fl2=glob(base+'map_ecc*'+survey+'_subchunk_2*.fits')
    m1=hp.ma(hp.read_map(fl1[0],[0,1,2]))
    m2=hp.ma(hp.read_map(fl2[0],[0,1,2]))
    tmask=wmask.astype(np.bool)|psmask|m1[1].mask|m2[1].mask
    for m in m1:
        m.mask=tmask
    for m in m2:
        m.mask=tmask
    cl=hp.anafast(m1,m2)
    l=np.arange(len(cl[0]))
    frac=(1-float(m1[0].mask.sum())/len(m1[0].mask))
    figure(4)
    plt.plot(l,1e12*cl[1]*l*(l+1)/(2*np.pi*frac),label='EE ' +survey)
    figure(5)
    plt.plot(l,1e12*cl[2]*l*(l+1)/(2*np.pi*frac),label='BB ' +survey)
    figure(6)
    plt.plot(l,1e12*cl[4]*l*(l+1)/(2*np.pi*frac),label='EB ' +survey)
       
figure(4)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EE halfRing cross-spectra for Elina ecc cal')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-20,20])
plt.savefig('cl_Elina ecc_cal_EE.png')

figure(5)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('BB halfRing cross-spectra for Elina ecc cal')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-20,20])
plt.savefig('cl_Elina ecc_cal_BB.png')
   
figure(6)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EB halfRing cross-spectra for Elina ecc cal')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-20,20])
plt.savefig('cl_Elina ecc_cal_EB.png')


figure(7)
figure(8)
figure(9)
quadse=['LFI18M','LFI19M','LFI20M']
quadsd=['18_23','19_22','20_21']

for quade,quadd in zip(quadse,quadsd):
    fl1eccs=glob(base+'map_eccs*'+quade+'*nominal*'+'_subchunk_1*.fits')
    fl2eccs=glob(base+'map_eccs*'+quade+'*nominal*'+'_subchunk_2*.fits')
    
    fl1ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_'+quadd+'_ringhalf_1_nominal_1s.fits']
    fl2ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_'+quadd+'_ringhalf_2_nominal_1s.fits']
    m1ecc=hp.ma(hp.read_map(fl1eccs[0],[0,1,2]))
    m2ecc=hp.ma(hp.read_map(fl2eccs[0],[0,1,2]))
    m1ddx9=hp.ma(hp.read_map(fl1ddx9[0],[0,1,2]))
    m2ddx9=hp.ma(hp.read_map(fl2ddx9[0],[0,1,2]))
    tmask=wmask.astype(np.bool)|psmask|m1ecc[1].mask|m2ecc[1].mask|m1ddx9[1].mask|m2ddx9[1].mask
    for m in m1ecc:
        m.mask=tmask
    for m in m2ecc:
        m.mask=tmask
    for m in m1ddx9:
        m.mask=tmask
    for m in m2ddx9:
        m.mask=tmask
        
    clecc=hp.anafast(m1ecc,m2ecc)
    clddx9=hp.anafast(m1ddx9,m2ddx9)
    l=np.arange(len(clecc[0]))
    fracecc=(1-float(m1ecc[0].mask.sum())/len(m1ecc[0].mask))
    fracddx9=(1-float(m1ddx9[0].mask.sum())/len(m1ddx9[0].mask))
    figure(7)
    hold(False)
    plt.plot(l,1e12*clecc[1]*l*(l+1)/(2*np.pi*fracecc),label='EE eccs')
    hold(True)
    plt.plot(l,1e12*clddx9[1]*l*(l+1)/(2*np.pi*fracddx9),label='EE ddx9')
    figure(8)
    hold(False)
    plt.plot(l,1e12*clecc[2]*l*(l+1)/(2*np.pi*fracecc),label='BB eccs')
    hold(True)
    plt.plot(l,1e12*clddx9[2]*l*(l+1)/(2*np.pi*fracddx9),label='BB ddx9')
    figure(9)
    hold(False)
    plt.plot(l,1e12*clecc[4]*l*(l+1)/(2*np.pi*fracecc),label='EB eccs')
    hold(True)
    plt.plot(l,1e12*clddx9[4]*l*(l+1)/(2*np.pi*fracddx9),label='EB ddx9')
           
    figure(7)
    xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EE halfRing cross-spectra for Elina ecc(samsmoothing) cal quad '+quadd)
    leg=plt.legend()
    plt.grid()
    plt.xlim([0,20])
    plt.ylim([-5,5])
    plt.savefig('cl_'+quadd+'_eccs_vs_ddx9_EE.png')

    figure(8)
    xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('BB halfRing cross-spectra for Elina ecc(samsmoothing) cal quad '+quadd)
    leg=plt.legend()
    plt.grid()
    plt.xlim([0,20])
    plt.ylim([-5,5])
    plt.savefig('cl_'+quadd+'_eccs_vs_ddx9_BB.png')

    figure(9)
    xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EB halfRing cross-spectra for Elina ecc(samsmoothing) cal quad '+quadd)
    leg=plt.legend()
    plt.grid()
    plt.xlim([0,20])
    plt.ylim([-5,5])
    plt.savefig('cl_'+quadd+'_eccs_vs_ddx9_EB.png')
    
figure(1)
figure(2)
figure(3)

fl1ecc200=glob(base+'map_ecc200*'+'*full*'+'_subchunk_1*.fits')
fl2ecc200=glob(base+'map_ecc200*'+'*full*'+'_subchunk_2*.fits')

fl1ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_ringhalf_1_full_1s.fits']
fl2ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120914_ringhalf_2_full_1s.fits']

m1ecc=hp.ma(hp.read_map(fl1ecc200[0],[0,1,2]))
m2ecc=hp.ma(hp.read_map(fl2ecc200[0],[0,1,2]))
m1ddx9=hp.ma(hp.read_map(fl1ddx9[0],[0,1,2]))
m2ddx9=hp.ma(hp.read_map(fl2ddx9[0],[0,1,2]))
tmask=wmask.astype(np.bool)|psmask|m1ecc[1].mask|m2ecc[1].mask|m1ddx9[1].mask|m2ddx9[1].mask
for m in m1ecc:
    m.mask=tmask
for m in m2ecc:
    m.mask=tmask
for m in m1ddx9:
    m.mask=tmask
for m in m2ddx9:
    m.mask=tmask
    
clecc=hp.anafast(m1ecc,m2ecc)
clddx9=hp.anafast(m1ddx9,m2ddx9)
l=np.arange(len(clecc[0]))
fracecc=(1-float(m1ecc[0].mask.sum())/len(m1ecc[0].mask))
fracddx9=(1-float(m1ddx9[0].mask.sum())/len(m1ddx9[0].mask))
figure(1)
hold(False)
plt.plot(l,1e12*clecc[1]*l*(l+1)/(2*np.pi*fracecc),label='EE ecc200')
hold(True)
plt.plot(l,1e12*clddx9[1]*l*(l+1)/(2*np.pi*fracddx9),label='EE ddx9')
figure(2)
hold(False)
plt.plot(l,1e12*clecc[2]*l*(l+1)/(2*np.pi*fracecc),label='BB ecc200')
hold(True)
plt.plot(l,1e12*clddx9[2]*l*(l+1)/(2*np.pi*fracddx9),label='BB ddx9')
figure(3)
hold(False)
plt.plot(l,1e12*clecc[4]*l*(l+1)/(2*np.pi*fracecc),label='EB ecc200')
hold(True)
plt.plot(l,1e12*clddx9[4]*l*(l+1)/(2*np.pi*fracddx9),label='EB ddx9')
       
figure(1)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EE halfRing cross-spectra for Elina ecc200 cal Full mission ')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_ecc200_vs_ddx9_EE.png')

figure(2)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('BB halfRing cross-spectra for Elina ecc200 cal Full mission ')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_ecc200_vs_ddx9_BB.png')

figure(3)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EB halfRing cross-spectra for Elina ecc200 cal Full mission')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_ecc200_vs_ddx9_EB.png')


figure(1)
figure(2)
figure(3)

fl1ecc200=glob(base+'map_ecc200*'+'*nominal*'+'_subchunk_1*.fits')
fl2ecc200=glob(base+'map_ecc200*'+'*nominal*'+'_subchunk_2*.fits')

fl1ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_1_nominal_1s.fits']
fl2ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_2_nominal_1s.fits']

m1ecc=hp.ma(hp.read_map(fl1ecc200[0],[0,1,2]))
m2ecc=hp.ma(hp.read_map(fl2ecc200[0],[0,1,2]))
m1ddx9=hp.ma(hp.read_map(fl1ddx9[0],[0,1,2]))
m2ddx9=hp.ma(hp.read_map(fl2ddx9[0],[0,1,2]))
tmask=wmask.astype(np.bool)|psmask|m1ecc[1].mask|m2ecc[1].mask|m1ddx9[1].mask|m2ddx9[1].mask
for m in m1ecc:
    m.mask=tmask
for m in m2ecc:
    m.mask=tmask
for m in m1ddx9:
    m.mask=tmask
for m in m2ddx9:
    m.mask=tmask
    
clecc=hp.anafast(m1ecc,m2ecc)
clddx9=hp.anafast(m1ddx9,m2ddx9)
l=np.arange(len(clecc[0]))
fracecc=(1-float(m1ecc[0].mask.sum())/len(m1ecc[0].mask))
fracddx9=(1-float(m1ddx9[0].mask.sum())/len(m1ddx9[0].mask))
figure(1)
hold(False)
plt.plot(l,1e12*clecc[1]*l*(l+1)/(2*np.pi*fracecc),label='EE ecc200')
hold(True)
plt.plot(l,1e12*clddx9[1]*l*(l+1)/(2*np.pi*fracddx9),label='EE ddx9')
figure(2)
hold(False)
plt.plot(l,1e12*clecc[2]*l*(l+1)/(2*np.pi*fracecc),label='BB ecc200')
hold(True)
plt.plot(l,1e12*clddx9[2]*l*(l+1)/(2*np.pi*fracddx9),label='BB ddx9')
figure(3)
hold(False)
plt.plot(l,1e12*clecc[4]*l*(l+1)/(2*np.pi*fracecc),label='EB ecc200')
hold(True)
plt.plot(l,1e12*clddx9[4]*l*(l+1)/(2*np.pi*fracddx9),label='EB ddx9')
       
figure(1)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EE halfRing cross-spectra for Elina ecc200 cal Nominal mission ')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_nominal_ecc200_vs_ddx9_EE.png')

figure(2)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('BB halfRing cross-spectra for Elina ecc200 cal Nominal mission ')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_nominal_ecc200_vs_ddx9_BB.png')

figure(3)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EB halfRing cross-spectra for Elina ecc200 cal Nominal mission')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_nominal_ecc200_vs_ddx9_EB.png')


figure(1)
figure(2)
figure(3)

fl1eccs=glob(base+'map_eccs*'+'*nominal*'+'_subchunk_1*.fits')
fl2eccs=glob(base+'map_eccs*'+'*nominal*'+'_subchunk_2*.fits')

fl1ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_1_nominal_1s.fits']
fl2ddx9=['/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/LFI_70_1024_20120912_ringhalf_2_nominal_1s.fits']

m1ecc=hp.ma(hp.read_map(fl1eccs[0],[0,1,2]))
m2ecc=hp.ma(hp.read_map(fl2eccs[0],[0,1,2]))
m1ddx9=hp.ma(hp.read_map(fl1ddx9[0],[0,1,2]))
m2ddx9=hp.ma(hp.read_map(fl2ddx9[0],[0,1,2]))
tmask=wmask.astype(np.bool)|psmask|m1ecc[1].mask|m2ecc[1].mask|m1ddx9[1].mask|m2ddx9[1].mask
for m in m1ecc:
    m.mask=tmask
for m in m2ecc:
    m.mask=tmask
for m in m1ddx9:
    m.mask=tmask
for m in m2ddx9:
    m.mask=tmask
    
clecc=hp.anafast(m1ecc,m2ecc)
clddx9=hp.anafast(m1ddx9,m2ddx9)
l=np.arange(len(clecc[0]))
fracecc=(1-float(m1ecc[0].mask.sum())/len(m1ecc[0].mask))
fracddx9=(1-float(m1ddx9[0].mask.sum())/len(m1ddx9[0].mask))
figure(1)
hold(False)
plt.plot(l,1e12*clecc[1]*l*(l+1)/(2*np.pi*fracecc),label='EE eccs')
hold(True)
plt.plot(l,1e12*clddx9[1]*l*(l+1)/(2*np.pi*fracddx9),label='EE ddx9')
figure(2)
hold(False)
plt.plot(l,1e12*clecc[2]*l*(l+1)/(2*np.pi*fracecc),label='BB eccs')
hold(True)
plt.plot(l,1e12*clddx9[2]*l*(l+1)/(2*np.pi*fracddx9),label='BB ddx9')
figure(3)
hold(False)
plt.plot(l,1e12*clecc[4]*l*(l+1)/(2*np.pi*fracecc),label='EB eccs')
hold(True)
plt.plot(l,1e12*clddx9[4]*l*(l+1)/(2*np.pi*fracddx9),label='EB ddx9')
       
figure(1)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EE halfRing cross-spectra for Elina ecc (sam smoothing) Nominal mission ')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_nominal_eccs_vs_ddx9_EE.png')

figure(2)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('BB halfRing cross-spectra for Elina ecc (sam smoothing) cal Nominal mission ')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_nominal_eccs_vs_ddx9_BB.png')

figure(3)
xlabel('Multipole'),ylabel('l(l+1)Cl, microK^2'),title('EB halfRing cross-spectra for Elina ecc (sam smoothing) cal Nominal mission')
leg=plt.legend()
plt.grid()
plt.xlim([0,20])
plt.ylim([-5,5])
plt.savefig('cl_70ghz_nominal_eccs_vs_ddx9_EB.png')






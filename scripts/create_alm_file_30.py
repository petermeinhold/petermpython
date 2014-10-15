import pyfits
import healpy as hp
from glob import glob
base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
basec='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/bandpass_correction/'

fl30nominal=glob(base+'/*30*nominal*.fits')
map30nomtqu=hp.read_map(fl30nominal[0],[0,1,2])
fl30nominalcorr=glob(basec+'/*correction*30*nominal*.fits')

corr30nomtqu=hp.read_map(fl30nominalcorr[0],[0,1,2])
corrmap30nom=[map30nomtqu[0]+corr30nomtqu[0],map30nomtqu[1]+corr30nomtqu[1],map30nomtqu[2]+corr30nomtqu[2]]
corrmap30nom_alm=hp.map2alm(corrmap30nom,pol=True)

hp.write_alm('alms/ddx9_30ghz_nom_alm_1024_t.fits',corrmap30nom_alm[0])
hp.write_alm('alms/ddx9_30ghz_nom_alm_1024_e.fits',corrmap30nom_alm[1])
hp.write_alm('alms/ddx9_30ghz_nom_alm_1024_b.fits',corrmap30nom_alm[2])

almt=pyfits.open('alms/ddx9_30ghz_nom_alm_1024_t.fits')
alme=pyfits.open('alms/ddx9_30ghz_nom_alm_1024_e.fits')
almb=pyfits.open('alms/ddx9_30ghz_nom_alm_1024_b.fits')

almt.append(alme[1])
almt.append(almb[1])
almt.writeto('alms/ddx9_30ghz_nom_alm_1024.fits')


fl70nominal=glob(base+'/*70*12_nominal*.fits')
map70nomtqu=hp.read_map(fl70nominal[0],[0,1,2])
fl70nominalcorr=glob(basec+'/*correction*70*nominal*.fits')

corr70nomtqu=hp.read_map(fl70nominalcorr[0],[0,1,2])
corrmap70nom=[map70nomtqu[0]+corr70nomtqu[0],map70nomtqu[1]+corr70nomtqu[1],map70nomtqu[2]+corr70nomtqu[2]]
corrmap70nom_alm=hp.map2alm(corrmap70nom,pol=True)

hp.write_alm('alms/ddx9_70ghz_nom_alm_1024_t.fits',corrmap70nom_alm[0])
hp.write_alm('alms/ddx9_70ghz_nom_alm_1024_e.fits',corrmap70nom_alm[1])
hp.write_alm('alms/ddx9_70ghz_nom_alm_1024_b.fits',corrmap70nom_alm[2])

almt=pyfits.open('alms/ddx9_70ghz_nom_alm_1024_t.fits')
alme=pyfits.open('alms/ddx9_70ghz_nom_alm_1024_e.fits')
almb=pyfits.open('alms/ddx9_70ghz_nom_alm_1024_b.fits')

almt.append(alme[1])
almt.append(almb[1])
almt.writeto('alms/ddx9_70ghz_nom_alm_1024.fits')



fl44nominal=glob(base+'/*44*_nominal*.fits')
map44nomtqu=hp.read_map(fl44nominal[0],[0,1,2])
fl44nominalcorr=glob(basec+'/*correction*44*nominal*.fits')

corr44nomtqu=hp.read_map(fl44nominalcorr[0],[0,1,2])
corrmap44nom=[map44nomtqu[0]+corr44nomtqu[0],map44nomtqu[1]+corr44nomtqu[1],map44nomtqu[2]+corr44nomtqu[2]]
corrmap44nom_alm=hp.map2alm(corrmap44nom,pol=True)

hp.write_alm('alms/ddx9_44ghz_nom_alm_1024_t.fits',corrmap44nom_alm[0])
hp.write_alm('alms/ddx9_44ghz_nom_alm_1024_e.fits',corrmap44nom_alm[1])
hp.write_alm('alms/ddx9_44ghz_nom_alm_1024_b.fits',corrmap44nom_alm[2])

almt=pyfits.open('alms/ddx9_44ghz_nom_alm_1024_t.fits')
alme=pyfits.open('alms/ddx9_44ghz_nom_alm_1024_e.fits')
almb=pyfits.open('alms/ddx9_44ghz_nom_alm_1024_b.fits')

almt.append(alme[1])
almt.append(almb[1])
almt.writeto('alms/ddx9_44ghz_nom_alm_1024.fits')

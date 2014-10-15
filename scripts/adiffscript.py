
from petermpython import planck_util as pu
import healpy as hp
from glob import glob

base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
fla=glob(base+'LFI*19_22*survey*.fits')
fla.sort()
flafull=glob(base+'LFI*19_22_full*.fits')

dbs1=pu.get_surv_diff_fname(fla[0],flafull[0],fwhm=1.,freq=70)
dbs2=pu.get_surv_diff_fname(fla[1],flafull[0],fwhm=1.,freq=70)
dbs3=pu.get_surv_diff_fname(fla[2],flafull[0],fwhm=1.,freq=70)
dbs4=pu.get_surv_diff_fname(fla[3],flafull[0],fwhm=1.,freq=70)
dbs5=pu.get_surv_diff_fname(fla[4],flafull[0],fwhm=1.,freq=70)

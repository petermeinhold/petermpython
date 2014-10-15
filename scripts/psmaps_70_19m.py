
from petermpython import planck_util as pu
import healpy as hp
from glob import glob
import quickring as qr
base='/global/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'
basechan='/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/Single_Radiometer/'
flachan19m=glob(basechan+'LFI*19M*survey*.fits')
flachan19m.sort()
flachan19mfull=glob(basechan+'LFI*19M_full*.fits')
dbs1=pu.get_surv_diff_fname(flachan19m[0],flachan19mfull[0],fwhm=1.,freq=70)
dbs2=pu.get_surv_diff_fname(flachan19m[1],flachan19mfull[0],fwhm=1.,freq=70)
dbs3=pu.get_surv_diff_fname(flachan19m[2],flachan19mfull[0],fwhm=1.,freq=70)
dbs4=pu.get_surv_diff_fname(flachan19m[3],flachan19mfull[0],fwhm=1.,freq=70)
dbs5=pu.get_surv_diff_fname(flachan19m[4],flachan19mfull[0],fwhm=1.,freq=70)

psdbs1=pu.rescan_to_psmap(dbs1,startring=3,stopring=5483,stepring=5,ang=qr.det_opening_angles['LFI19M'])
psdbs2=pu.rescan_to_psmap(dbs2,startring=5484,stopring=10957,stepring=5,ang=qr.det_opening_angles['LFI19M'])
psdbs3=pu.rescan_to_psmap(dbs3,startring=10958,stopring=16454,stepring=5,ang=qr.det_opening_angles['LFI19M'])
psdbs4=pu.rescan_to_psmap(dbs4,startring=16455,stopring=21482,stepring=5,ang=qr.det_opening_angles['LFI19M'])
psdbs5=pu.rescan_to_psmap(dbs5,startring=21483,stopring=27404,stepring=5,ang=qr.det_opening_angles['LFI19M'])
imshow(psdbs1)
title('Pseudomap, LFI 70 GHz 19M SS1-full')
xlabel('Clock angle, 0-250/rev'),ylabel('relative PID/5')
savefig('plots/oct10/psmap19m_s1_full.png')

imshow(psdbs2)
title('Pseudomap, LFI 70 GHz 19M SS2-full')
xlabel('Clock angle, 0-250/rev'),ylabel('relative PID/5')
savefig('plots/oct10/psmap19m_s2_full.png')

imshow(psdbs3)
title('Pseudomap, LFI 70 GHz 19M SS3-full')
xlabel('Clock angle, 0-250/rev'),ylabel('relative PID/5')
savefig('plots/oct10/psmap19m_s3_full.png')

imshow(psdbs4)
title('Pseudomap, LFI 70 GHz 19M SS4-full')
xlabel('Clock angle, 0-250/rev'),ylabel('relative PID/5')
savefig('plots/oct10/psmap19m_s4_full.png')

imshow(psdbs5)
title('Pseudomap, LFI 70 GHz 19M SS5-full')
xlabel('Clock angle, 0-250/rev'),ylabel('relative PID/5')
savefig('plots/oct10/psmap19m_s5_full.png')

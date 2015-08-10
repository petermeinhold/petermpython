import sys
sys.path.append('/global/homes/p/peterm/petermpython')
sys.path.append('/global/homes/p/peterm/paperplots/python/scripts')
import os
import numpy as np
import time
from planck import Planck
from planck.metadata import latest_exchange
from glob import glob
import pyfits
import testenv
from planck.metadata import latest_exchange
import datetime
from planck_util_prm import get_ucds

now=datetime.datetime.now()
formatted_date= now.strftime("%Y%m%d")
horn=int(sys.argv[1])
freq=70
if horn>23:
	freq=44
if horn>26:
	freq=30
chans=['LFI%sM' %str(horn),'LFI%sS' %str(horn)]

gmf={}
print('Frequency',freq)	
ucds=get_ucds(db='/global/homes/p/peterm/ucds-dx11-delta.db', freq=freq)
keys=ucds.keys()
for key in keys:
	for d in ['0','1']:
		uchan=key+d
		print(uchan)
		dchan=uchan.replace('S','S1')
		dchan=dchan.replace('M','M0')
		gmf[dchan]=np.mean(ucds[key]['vsky'+d])/np.mean(ucds[key]['vref'+d])

conf={}
conf['cal_tag']='usp_1'
pl=Planck()

weights=np.genfromtxt('/global/homes/p/peterm/weights_march15_2010.txt',dtype=str,skip_header=1,autostrip=True,delimiter=',')
diodeweights={}
for weight in weights:
	diodeweights[weight[0]]=float(weight[1])

def get_raw_diff(chan,targetod,effobt,gmf,diodeweights,rawdir='/global/scratch2/sd/planck/user/peterm/rawdata/spike_corrected/'):
	#function to grab 3 OD's of raw data (+-1 od), cut to match eff obt, difference and return
	horn=chan[3:5]
	rad=chan[-1]
	diodes={}
	diodes['M']=['00','01']
	diodes['S']=['10','11']
	diff=np.zeros(len(effobt))
	for dnum,diode in enumerate(diodes[rad]):
		sky=np.array([])
		ref=np.array([])
		rawobt=np.array([])
		for od in [targetod-1,targetod,targetod+1]:
			rawfile=rawdir+'/od%s_rca%s%s.fits' %(str(od),horn,diode)
			if not(os.path.exists(rawfile)):
				print('no raw file',rawfile)
				stat=False
				return stat,diff
			hdulist=pyfits.open(rawfile)
			rawobt=np.concatenate([rawobt,hdulist[1].data['TIME']])
			sky=np.concatenate([sky,hdulist[1].data['SKY']])
			ref=np.concatenate([ref,hdulist[1].data['REF']])
		sky=sky[(rawobt>=effobt[0]) & (rawobt<=effobt[-1])]
		ref=ref[(rawobt>=effobt[0]) & (rawobt<=effobt[-1])]
		#cut to match effobt
		if dnum==0:
			rawobt=rawobt[(rawobt>=effobt[0]) & (rawobt<=effobt[-1])]
			rawdiff = np.interp(effobt,rawobt,diodeweights[chan+'-'+diode]*(sky-ref*gmf[chan+diode]))
		else:
			rawobt=rawobt[(rawobt>=effobt[0]) & (rawobt<=effobt[-1])]
			rawdiff += np.interp(effobt,rawobt,diodeweights[chan+'-'+diode]*(sky-ref*gmf[chan+diode]))
			diff=rawdiff
		stat=True
	if np.isnan(diff).sum() != 0:
		stat=False    
	return stat,diff
	
badlist=[]

for od in range(91,1605):
	print(od)
	if od==1540:
		continue
	filename = '%s%03d-%04d-%s.fits' % ("L", freq, od, "C")
	output_folder = "/global/scratch2/sd/planck/user/peterm/data/%s/%04d/" % (conf["cal_tag"], od)
	if os.path.isfile(output_folder+filename):
		continue
	efftype = "C"
	eff_infilename = latest_exchange(freq, ods=od, exchangefolder = "/project/projectdirs/planck/data/mission/lfi_ops_dx11_delta/", type = efftype)
	with pyfits.open(eff_infilename) as eff_infile:
		effobt=eff_infile['OBT'].data['OBT']
		for chan in chans:
			stat,diffdata=get_raw_diff(chan,od,effobt,gmf,diodeweights)
			if stat==False:
				print('some Nans, od= ',od)
			eff_infile[chan].data[chan] = diffdata 
		if stat==False:
			badlist.append(od)
			continue
		try:
			os.mkdir(output_folder)
		except:
			pass
		eff_infile.writeto(output_folder + filename, clobber=True)
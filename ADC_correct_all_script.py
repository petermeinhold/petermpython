from numpy import *
import pyfits as fits
import cPickle as pic
import scipy.interpolate.fitpack as fit
import sys
import shutil
import os
from glob import glob

adc_dir='/global/scratch2/sd/planck/user/peterm/DPC_ADC_DX12/'
rawdata_dir='/global/scratch2/sd/planck/user/peterm/rawdata/DDS/fits_toi/'
adc_corr_data_dir='/global/scratch2/sd/planck/user/peterm/rawdata/adc_12_corrected/'


class LFI_response:
	""" Class to store response curve in similar format to
	LFI_response data object
	"""
	def __init__(self,keys,sky_Vi,sky_Vo,ref_Vi,ref_Vo):
		self.keys = keys             # Dictionary of keys and values
		self.sky_volt_in = sky_Vi
		self.sky_volt_out = sky_Vo
		self.load_volt_in = ref_Vi
		self.load_volt_out = ref_Vo
	  
def read_LFI_response(filen):
	pf = open(filen,'rb')
	resp = pic.load(pf)
	pf.close()
	return resp

horn=sys.argv[1]

diodes=[]
for rad in ['0','1']:
	for diode in ['0','1']:
		diodes.append(str(horn)+rad+diode)

for diode in diodes:
	#first do 91-952, read response file then go thru ODs
	fl=glob(adc_dir+'/'+'ADCDX12_%s*_952_*.pic' % diode)
	if len(fl)>0:
		respfile=fl[0]
		if os.path.exists(respfile):
			resp = read_LFI_response(respfile)
			horn = resp.keys['horn']
			rad  = resp.keys['radiometer']
			det  = resp.keys['detector'] 
			sky_spline = fit.splrep(resp.sky_volt_out,resp.sky_volt_in,s=0.0)
			ref_spline = fit.splrep(resp.load_volt_out,resp.load_volt_in,s=0.0)
			for od in range(91,953):
				rawfile=rawdata_dir+'od%s_rca%s.fits' % (od,diode)
				corrfile=adc_corr_data_dir+'od%s_rca%s.fits' % (od,diode)
				if not(os.path.exists(corrfile)):
					if os.path.exists(rawfile):
						hdulist = fits.open(rawfile)
						data = hdulist[1].data
						sky_volt  = data.field('SKY')
						load_volt = data.field('REF')						
						sky_cor = fit.splev(sky_volt,sky_spline)
						load_cor = fit.splev(load_volt,ref_spline)
						data.field('SKY')[:] = sky_cor 
						data.field('REF')[:] = load_cor
						hdulist[1].header.update('hierarch instrument','LFI_ADC_corr')
						hdulist[1].header.update('hierarch ADC_correction','True')
						print 'ADC correcting %s -> %s' % (rawfile,corrfile)
						hdulist.writeto(corrfile)
	if len(fl)==0:
		for od in range(91,953):
			rawfile=rawdata_dir+'od%s_rca%s.fits' % (od,diode)
			corrfile=adc_corr_data_dir+'od%s_rca%s.fits' % (od,diode)
			if os.path.exists(rawfile):
				if not(os.path.exists(corrfile)):
					shutil.copy(rawfile,corrfile)
		
		#now do the od's later than 952 - in the case there is only 1 response file, use it
	fl=glob(adc_dir+'/'+'ADCDX12_%s*_953_*.pic' % diode)
	if len(fl)>0:
		respfile=fl[0]
		if os.path.exists(respfile):
			resp = read_LFI_response(respfile)
			horn = resp.keys['horn']
			rad  = resp.keys['radiometer']
			det  = resp.keys['detector'] 
			sky_spline = fit.splrep(resp.sky_volt_out,resp.sky_volt_in,s=0.0)
			ref_spline = fit.splrep(resp.load_volt_out,resp.load_volt_in,s=0.0)
			for od in range(952,1605):
				rawfile=rawdata_dir+'od%s_rca%s.fits' % (od,diode)
				corrfile=adc_corr_data_dir+'od%s_rca%s.fits' % (od,diode)
				if not(os.path.exists(corrfile)):
					if os.path.exists(rawfile):
						hdulist = fits.open(rawfile)
						data = hdulist[1].data
						sky_volt  = data.field('SKY')
						load_volt = data.field('REF')						
						sky_cor = fit.splev(sky_volt,sky_spline)
						load_cor = fit.splev(load_volt,ref_spline)
						data.field('SKY')[:] = sky_cor 
						data.field('REF')[:] = load_cor
						hdulist[1].header.update('hierarch instrument','LFI_ADC_corr')
						hdulist[1].header.update('hierarch ADC_correction','True')
						print 'ADC correcting %s -> %s' % (rawfile,corrfile)
						hdulist.writeto(corrfile)
	if len(fl)==0:
		for od in range(952,1605):
			rawfile=rawdata_dir+'od%s_rca%s.fits' % (od,diode)
			corrfile=adc_corr_data_dir+'od%s_rca%s.fits' % (od,diode)
			if os.path.exists(rawfile):
				if not(os.path.exists(corrfile)):
					shutil.copy(rawfile,corrfile)
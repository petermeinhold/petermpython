import numpy as np
import pyfits as fits
import sys
import os

adc_dir='/global/scratch2/sd/planck/user/peterm/ADC_correction/'
rawdata_dir='/global/scratch2/sd/planck/user/peterm/rawdata/DDS/fits_toi/'
adc_corr_data_dir='/global/scratch2/sd/planck/user/peterm/rawdata/adc_12_corrected/'
spike_corr_data_dir='/global/scratch2/sd/planck/user/peterm/rawdata/spike_corrected/'


def get_spike_template(obt,sky,ref):
    """
    function finds unique fractional second bins for a bolus of obt,sky,ref data 
    bins sky and ref to those bins, subtracts offsets, returns sbins and templates
    """
    obtt=obt/65536.   #convert the OBT to seconds
    obttf,obtti=np.modf(obtt)   #find the fractional part, so we can bin to a consistent phase
    sbins=np.unique(obttf)   #find unique values of the fractional time for this day and chan, will be the template bins
    #do the binning
    rtemplate=np.zeros(len(sbins))
    stemplate=np.zeros(len(sbins))
    for t,sbin in enumerate(sbins):
        rtemplate[t]=np.mean(ref[obttf==sbin])
        stemplate[t]=np.mean(sky[obttf==sbin])
    return stemplate-np.mean(stemplate),rtemplate-np.mean(rtemplate) 
    
def subtract_spike_template(obt,sky,ref,stemplate,rtemplate):
    """
    function finds unique fractional second bins for a bolus of obt,sky,ref data 
    bins sky and ref to those bins, subtracts offsets, returns sbins and templates
    """
    obtt=obt/65536.   #convert the OBT to seconds
    obttf,obtti=np.modf(obtt)   #find the fractional part, so we can bin to a consistent phase
    sbins=np.unique(obttf)   #find unique values of the fractional time for this day and chan, will be the template bins
    subref=ref.copy()
    subsky=sky.copy()
    #do the subtraction
    for t,sbin in enumerate(sbins):
        subsky[obttf==sbin] = sky[obttf==sbin]-stemplate[t]
        subref[obttf==sbin] = ref[obttf==sbin]-rtemplate[t]
    return subsky,subref
    
horn=sys.argv[1]
diodes=[]
for rad in ['0','1']:
	for diode in ['0','1']:
	    diodes.append(str(horn)+rad+diode)
print(diodes)

for diode in diodes:
    for od in range(91,1603):
        rawfile=adc_corr_data_dir+'od%s_rca%s.fits' % (od,diode)
        corrfile=spike_corr_data_dir+'od%s_rca%s.fits' % (od,diode)
        if not(os.path.exists(corrfile)):
            if os.path.exists(rawfile):
                hdulist = fits.open(rawfile)
                data = hdulist[1].data
                sky  = data.field('SKY')
                ref = data.field('REF')
                obt=data.field('OBT')
                sky_template,ref_template=get_spike_template(obt,sky,ref)	
                sub_sky,sub_ref=subtract_spike_template(obt,sky,ref,sky_template,ref_template)
                data.field('SKY')[:] = sub_sky 
                data.field('REF')[:] = sub_ref
                hdulist[1].header.update('hierarch instrument','LFI_spikes_removed')
                hdulist[1].header.update('hierarch spikes_removed','True')
                print('1 Hz Spike removal %s -> %s' % (rawfile,corrfile)
                hdulist.writeto(corrfile)

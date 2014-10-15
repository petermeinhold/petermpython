import numpy as np
from planck import private
import h5py
import sqlite3
import scikits.statsmodels.api as sm
from planck.Planck import Planck
import planck_util as pu

MS = 'MS'

#conn = sqlite3.connect("ucdb_dx9_30.db")
conn = sqlite3.connect("/project/projectdirs/planck/user/zonca/issues/deltav/ucds-dx9.db")
c = conn.cursor()
ch_data_dtype = np.dtype(
        [('pID',np.long), ('vref',np.double),('vref0',np.double),('vref1',np.double), 
        ('vsky0',np.double),('vsky1',np.double), 
        ('gain',np.double), ('nominal_gain',np.double), ('dipoleT',np.double), 
        ('l14k',np.double), ('ts1l',np.double), ('ts2l',np.double), ('ts3l',np.double), 
        ('ts4l',np.double), ('ts5l',np.double), ('ts6l',np.double), ('ts1r',np.double), 
        ('ts2r',np.double), ('ts3r',np.double), ('ts4r',np.double), ('ts5r',np.double), 
        ('ts6r',np.double), ('beusvc',np.double), ('beusci1',np.double), ('beusci2',np.double), 
        ('beusci3',np.double), ('beusci4',np.double), ('lbem1',np.double), ('lbem2',np.double), 
        ('lfem1',np.double), ('lfem2',np.double), ('ldaq1',np.double), ('rbem1',np.double), 
        ('rbem2',np.double), ('rfem1',np.double), ('rfem2',np.double), ('rdaq1',np.double)])

#def extract_fields(structarray, fields):
#    """Extract fields as a matrix"""
#    return structarray[fields].view(np.double).reshape((-1, len(fields)))
#
#def dvv_fit(gains, dvref, weights, add_constant=False):
#    if add_constant:
#        model = sm.WLS(gains, sm.add_constant(dvref), weights=weights)
#    else:
#        model = sm.WLS(gains, (dvref), weights=weights)
#    return model.fit().params

#def write_gain(caltype, freq, cal):
#    h5filename = "~/cals/dstcal/cal/%s_%d.h5" % (caltype, freq)
#    with h5py.File(h5filename, mode='w') as f:
#    for i,ch in enumerate(cal.dtype.names):
#       flat_data[:, i] = cal[ch]
#    with h5py.File(h5filename.replace('.h5','_flat.h5'), mode='w') as h5f:
#        h5f.create_dataset('data', data=flat_data.flatten())

#base of deltav is dpc-style cal with full mission, fit first iteration and temperature destriping
#dipole_cal = h5py.File("/project/projectdirs/planck/user/zonca/issues/dstcal/DX8S_dm_ff_td_smt.h5", mode='r')['data'][:]
#def zero_cal():
#    return np.zeros(len(dipole_cal), dtype=dipole_cal.dtype) 
#dvv_cal = zero_cal()
#allpids=h5py.File('all_pids_030.h5')['data'][:len(dipole_cal)]

freq = 30
f = Planck().f[freq]

#construct set of hskp data vectors to fit to:
all_data = {}
for ch in f.ch:
    print(ch)
    horn = ch.horn
    rad = ch.n
    chtag = "LFI%d%s" % (horn, MS[rad])
    c.execute(' '.join([ "select sci.pointingID as pID,",
             "sci.mean_ref as vref,",
             "sci_d0.mean_ref as vref0,",
             "sci_d1.mean_ref as vref1,",
             "sci_d0.mean_sky as vsky0,",
             "sci_d1.mean_sky as vsky1,",
             "rg.gain as gain,",
             "sg.gain as nominal_gain,",
             "dipole.deltaT as dipoleT,",
             "hfihk.l1_a5 as l14k,",
             "lfihk.feu_left_side_wall as ts1l,",
             "lfihk.feu_cold_plate_left_inner as ts2l,",
             "lfihk.feu_cold_plate_right_inner as ts3l,",
             "lfihk.feu_left_bottom_fh25 as ts4l,",
             "lfihk.feu_cold_plate_far_left as ts5l,",
             "lfihk.feu_cone_left_part as ts6l,",
             "lfihk.feu_right_bottom_fh28 as ts1r,",
             "lfihk.feu_cone_right_part as ts2r,",
             "lfihk.feu_right_side_wall as ts3r,",
             "lfihk.feu_cold_plate_far_right as ts4r,",
             "lfihk.feu_fh28_flange as ts5r,",
             "lfihk.feu_right_bottom_fh26 as ts6r,",
             "lfihk.beu_service as beusvc,",
             "lfihk.beu_science1 as beusci1,",
             "lfihk.beu_science2 as beusci2,",
             "lfihk.beu_science3 as beusci3,",
             "lfihk.beu_science4 as beusci4,",
             "lfihk.beu_l_bem1 as lbem1,",
             "lfihk.beu_l_bem2 as lbem2,",
             "lfihk.beu_l_fem1 as lfem1,",
             "lfihk.beu_l_fem2 as lfem2,",
             "lfihk.beu_r_bem1 as rbem1,",
             "lfihk.beu_r_bem2 as rbem2,",
             "lfihk.beu_r_fem1 as rfem1,",
             "lfihk.beu_r_fem2 as rfem2,",
             "lfihk.beu_l_daq1 as ldaq1,",
             "lfihk.beu_r_daq1 as rdaq1",
             "from sci%d%d_weighted as sci" % (horn, rad),
             "join sci%d%d%d as sci_d0 using(pointingID)" % (horn,rad,0),
             "join sci%d%d%d as sci_d1 using(pointingID)" % (horn,rad,1),
             "join raw_gains%d%d as rg using(pointingID)"% (horn, rad),
             "join smoothed_gains%d%d as sg using(pointingID)" % (horn, rad),
             "join dipole_%d%d as dipole using(pointingID)" % (horn, rad),
             "join hfi_hk_temperatures as hfihk using(pointingID)",
             "join lfi_hk_temperatures as lfihk using(pointingID)",
             "where vref > 0",
             "and rg.gain > 0",
             "and (sci.pointingID < 10911 or sci.pointingID > 10920)",
             "order by sci.pointingID" ]))
    all_list=[list(a) for a in c]
    for a in all_list:
        for i in range(len(a)):
            if a[i]=='':
                a[i]=0
    
    ch_data = np.array([tuple(a) for a 
    in all_list],dtype=ch_data_dtype)
    all_data[ch.tag] = ch_data


fitlist=['lbem1', 'ts1r']
fitlist=['ts1l','ts1r','ts4r','ts6r','beusvc','lbem1','ldaq1']
ranges = [slice(0, 21522)]

ranges = []
for surv in [1,2,3,4]:
    ranges.append(slice(private.survey[surv].SPID[0], min(private.survey[surv].SPID[1]+1,len(allpids)-1)))

#ranges = [slice(0, 10988), slice(10988, 21522)] 
ranges = [slice(0, 21522)] 

for rang in ranges:
    dpc_rang = slice(all_data[ch.tag]['pID'].searchsorted(allpids['pointID_unique'][rang.start]), ch_data['pID'].searchsorted(allpids['pointID_unique'][rang.stop]))
    for ch in f.ch:
        print(ch)
        hkvecs = extract_fields(all_data[ch.tag], fitlist)
        fit_result = dvv_fit(all_data[ch.tag]['gain'][dpc_rang], hkvecs[dpc_rang], all_data[ch.tag]['dipoleT'][dpc_rang], add_constant=True)
        print(fit_result)
        hkmodel=np.sum([fit_result[i]*hkvecs[:,i] for  i in range(len(fitlist))],axis=0) + fit_result[-1]
        dvv_cal[ch.tag][rang] = np.interp(allpids['pointID_unique'][rang], ch_data['pID'][dpc_rang], hkmodel[dpc_rang])
        
write_gain("hk7", freq, dvv_cal)

import planck_util_prm as pu
from glob import glob
from scipy.stats import nanmean,nanstd,nanmedian
ucds30_ddx9=pu.get_ucds(freq=30)
ucds44_ddx9=pu.get_ucds(freq=44)
ucds70_ddx9=pu.get_ucds(freq=70)

ucds30_dx10=pu.get_ucds(db="/global/homes/p/peterm/ucds-dx10.db",freq=30)
ucds44_dx10=pu.get_ucds(db="/global/homes/p/peterm/ucds-dx10.db",freq=44)
ucds70_dx10=pu.get_ucds(db="/global/homes/p/peterm/ucds-dx10.db",freq=70)

ucds_ddx9=dict(ucds70_ddx9,**ucds44_ddx9)
ucds_ddx9.update(ucds30_ddx9)
ucds_dx10=dict(ucds70_dx10,**ucds44_dx10)
ucds_dx10.update(ucds30_dx10)

rcalist=[]
for r in range(18,29):
    for d in ['M','S']:
        rcalist.append('LFI'+str(r)+d)
    

for rca in rcalist:
    figure()
    dvref0_ddx9=1/(ucds_ddx9[rca]['vsky0']-ucds_ddx9[rca]['vref0'])
    dvref1_ddx9=1/(ucds_ddx9[rca]['vsky1']-ucds_ddx9[rca]['vref1'])
    plot(ucds_ddx9[rca]['pID'],dvref0_ddx9/np.mean(dvref0_ddx9),label='vsky-vref diode0 dDX9')
    plot(ucds_ddx9[rca]['pID'],dvref1_ddx9/np.mean(dvref1_ddx9),label='vsky-vref diode1 dDX9')

    dvref0_dx10=1/(ucds_dx10[rca]['vsky0']-ucds_dx10[rca]['vref0'])
    dvref1_dx10=1/(ucds_dx10[rca]['vsky1']-ucds_dx10[rca]['vref1'])
    plot(ucds_dx10[rca]['pID'],dvref0_dx10/np.mean(dvref0_dx10),label='vsky-vref diode0 dx10')
    plot(ucds_dx10[rca]['pID'],dvref1_dx10/np.mean(dvref1_dx10),label='vsky-vref diode1 dx10')

    title('Normalized diode inverse Sky-Ref, ' + rca)
    xlabel('pID')
    ylabel('Relative value')
    leg=legend()
    pu.thicklegendlines(leg)
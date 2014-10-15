cal1='dx9'
base='/global/project/projectdirs/planck/data/mission/DPC_maps/'+cal1+'/lfi/'
fl=glob(base+'*70*7_survey*.fits')
fl
qdiff12=pu.get_surv_diff_fname(fl[0],fl[1],pol='Q',freq=70)
qdiff23=pu.get_surv_diff_fname(fl[1],fl[2],pol='Q',freq=70)
qdiff13=pu.get_surv_diff_fname(fl[0],fl[2],pol='Q',freq=70)
psm12s2=pu.rescan_to_psmap(qdiff12,startring=5484,stopring=10957)
psm23s2=pu.rescan_to_psmap(qdiff23,startring=5484,stopring=10957)
psm12s1=pu.rescan_to_psmap(qdiff12,startring=3,stopring=5483)
psm23s1=pu.rescan_to_psmap(qdiff23,startring=3,stopring=5483)
psm13s1=pu.rescan_to_psmap(qdiff13,startring=3,stopring=5483)

figure('DX9 70 GHz SS1-SS3 rescan SS1')
pcolor transpose(psm13s1),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz SS1-SS2 rescan SS1')
pcolor transpose(psm12s1),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz SS1-SS2 rescan SS2')
pcolor transpose(psm12s2),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz SS2-SS3 rescan SS1')
pcolor transpose(psm23s1),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz SS2-SS3 rescan SS2')
pcolor transpose(psm23s2),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')

fl=glob(base+'*70*20_*_survey*.fits')
q2021diff12=pu.get_surv_diff_fname(fl[0],fl[1],pol='Q',freq=70)
q2021diff23=pu.get_surv_diff_fname(fl[1],fl[2],pol='Q',freq=70)
q2021diff13=pu.get_surv_diff_fname(fl[0],fl[2],pol='Q',freq=70)
ps2021m12s2=pu.rescan_to_psmap(q2021diff12,startring=5484,stopring=10957)
ps2021m23s2=pu.rescan_to_psmap(q2021diff23,startring=5484,stopring=10957)
ps2021m12s1=pu.rescan_to_psmap(q2021diff12,startring=3,stopring=5483)
ps2021m23s1=pu.rescan_to_psmap(q2021diff23,startring=3,stopring=5483)
ps2021m13s1=pu.rescan_to_psmap(q2021diff13,startring=3,stopring=5483)
ps2021m13s3=pu.rescan_to_psmap(q2021diff13,startring=10958,stopring=16454)

figure('DX9 70 GHz Q 2021 SS1-SS3 rescan SS3')
pcolor transpose(ps2021m13s3),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')

figure('DX9 70 GHz Q 2021 SS1-SS3 rescan SS1')
pcolor transpose(ps2021m13s1),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz Q 2021 SS1-SS2 rescan SS1')
pcolor transpose(ps2021m12s1),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz Q 2021 SS1-SS2 rescan SS2')
pcolor transpose(ps2021m12s2),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz Q 2021 SS2-SS3 rescan SS1')
pcolor transpose(ps2021m23s1),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')
figure('DX9 70 GHz Q 2021 SS2-SS3 rescan SS2')
pcolor transpose(ps2021m23s2),vmin=-10e-6,vmax=10e-6,cmap=get_cmap('gray')




dvcals=glob('cals/*dvv*.txt')
dx8cals=glob('cals/*dx8*.txt')
dvv18=pu.import_csv_cals(dvcals[0])
dvv20=pu.import_csv_cals(dvcals[2])
dvv21=pu.import_csv_cals(dvcals[3])
dvv22=pu.import_csv_cals(dvcals[4])
dvv23=pu.import_csv_cals(dvcals[5])

figure()
plot dvv19[0,:],dvv19[1,:]/mean(dvv19[1,:]),label='19 SS1'
plot dvv19[0,:]-10957,dvv19[1,:]/mean(dvv19[1,:]),label='19 SS3'
plot dvv20[0,:],dvv20[1,:]/mean(dvv20[1,:]),label='20 SS1'
plot dvv20[0,:]-10957,dvv20[1,:]/mean(dvv20[1,:]),label='20 SS3'
plot dvv21[0,:],dvv21[1,:]/mean(dvv21[1,:]),label='21 SS1'
plot dvv21[0,:]-10957,dvv21[1,:]/mean(dvv21[1,:]),label='21 SS3'
plot dvv22[0,:],dvv22[1,:]/mean(dvv22[1,:]),label='22 SS1'
plot dvv22[0,:]-10957,dvv22[1,:]/mean(dvv22[1,:]),label='22 SS3'
plot dvv23[0,:],dvv23[1,:]/mean(dvv23[1,:]),label='23 SS1'
plot dvv23[0,:]-10957,dvv23[1,:]/mean(dvv23[1,:]),label='23 SS3'
leg=legend()
pu.thicklegendlines(leg)

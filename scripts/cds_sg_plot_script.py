import planck_util_prm as pu

t4k=pu.get_cds_key_value()
tfh28=pu.get_cds_key_value(key='feu_temp',value='fh28_flange')
tlbem1=pu.get_cds_key_value(key='beu_temp',value='l_bem1')

sgt4k=pu.spectrogram_cds(t4k)
sgfh28=pu.spectrogram_cds(tfh28)
sglbem1=pu.spectrogram_cds(tlbem1)

figure(figsize=[10,4])
pcolormesh(np.array(sgt4k['od']).T,np.array(sgt4k['frequency']).T,log(sgt4k['spectrogram']).T,vmax=3,vmin=-14)
title('Spectrogram, 4 K stage temperature (cernox)')
xlabel('Operational Day'),ylabel('Frequency, mHz')
xlim([91,1548]),ylim([0,8.3])
cb=colorbar()
cb.set_label('$log(K^2 / Hz)$')
show()
savefig('/global/homes/p/peterm/plots/cds/cds_sg_t4k.png')

figure(figsize=[10,4])
pcolormesh(np.array(sglbem1['od']).T,np.array(sglbem1['frequency']).T,log(sglbem1['spectrogram']).T,vmax=3,vmin=-14)
title('Spectrogram, Backend Unit temperature (l_bem1)')
xlabel('Operational Day'),ylabel('Frequency, mHz')
xlim([91,1548]),ylim([0,8.3])
cb=colorbar()
cb.set_label('$log(K^2 / Hz)$')
show()
savefig('/global/homes/p/peterm/plots/cds/cds_sg_beu.png')

figure(figsize=[10,4])
pcolormesh(np.array(sgfh28['od']).T,np.array(sgfh28['frequency']).T,log(sgfh28['spectrogram']).T,vmax=3,vmin=-14)
title('Spectrogram, Front End Unit temperature (feedhorn 28 flange)')
xlabel('Operational Day'),ylabel('Frequency, mHz')
xlim([91,1548]),ylim([0,8.3])
cb=colorbar()
cb.set_label('$log(K^2 / Hz)$')
show()
savefig('/global/homes/p/peterm/plots/cds/cds_sg_feu.png')

cds=pu.get_cds()

chans=[]
for h in range(18,29):
    for ext in ['00','01','10','11']:
        chans.append('rca'+np.str(h)+ext)
        
for chan in chans[36:]:
    print chan
    sg=pu.spectrogram_cds(cds[chan]['average_dif'])
    figure(figsize=[10,4])
    if np.int(chan[3:5])<24:
        pcolormesh(np.array(sg['od']).T,np.array(sg['frequency']).T,log(sg['spectrogram']).T,vmax=-14,vmin=-20)
    if (np.int(chan[3:5])<27) and (np.int(chan[3:5])>23):
        pcolormesh(np.array(sg['od']).T,np.array(sg['frequency']).T,log(sg['spectrogram']).T,vmax=-18,vmin=-24)        
    if np.int(chan[3:5])>26:
        pcolormesh(np.array(sg['od']).T,np.array(sg['frequency']).T,log(sg['spectrogram']).T,vmax=-14,vmin=-20)
    title('Spectrogram, average sky-ref '+chan)
    xlabel('Operational Day'),ylabel('Frequency, mHz')
    xlim([91,1548]),ylim([0,8.3])
    cb=colorbar()
    cb.set_label('$log(V^2 / Hz)$')
    show()
    savefig('/global/homes/p/peterm/plots/cds/cds_sg_'+chan+'.png')




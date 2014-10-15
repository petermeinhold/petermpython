
m70y1ebsm10deg=pu.get_dx10_map(eb=True,surv='yr1',fwhm=np.radians(10.0))
m70y2ebsm10deg=pu.get_dx10_map(eb=True,surv='yr2',fwhm=np.radians(10.0))
m70y3ebsm10deg=pu.get_dx10_map(eb=True,surv='yr3',fwhm=np.radians(10.0))


hp.mollview m70y1ebsm10deg[1],min=-2e-6,max=2e-6,title='DX10 Year 1, 70 GHz E map, smoothed 10 degrees',unit='K'
savefig('emap_yr1_70ghz_10deg.png')

hp.mollview m70y1ebsm10deg[2],min=-2e-6,max=2e-6,title='DX10 Year 1, 70 GHz B map, smoothed 10 degrees',unit='K'
savefig('bmap_yr1_70ghz_10deg.png')

hp.mollview m70y2ebsm10deg[1],min=-2e-6,max=2e-6,title='DX10 Year 2, 70 GHz E map, smoothed 10 degrees',unit='K'
savefig('emap_yr2_70ghz_10deg.png')

hp.mollview m70y2ebsm10deg[2],min=-2e-6,max=2e-6,title='DX10 Year 2, 70 GHz B map, smoothed 10 degrees',unit='K'
savefig('bmap_yr2_70ghz_10deg.png')

hp.mollview m70y3ebsm10deg[1],min=-2e-6,max=2e-6,title='DX10 Year 3, 70 GHz E map, smoothed 10 degrees',unit='K'
savefig('emap_yr3_70ghz_10deg.png')

hp.mollview m70y3ebsm10deg[2],min=-2e-6,max=2e-6,title='DX10 Year 3, 70 GHz B map, smoothed 10 degrees',unit='K'
savefig('bmap_yr3_70ghz_10deg.png')


m70y1qusm10deg=pu.get_dx10_map(surv='yr1',fwhm=np.radians(10.0))
m70y2qusm10deg=pu.get_dx10_map(surv='yr2',fwhm=np.radians(10.0))
m70y3qusm10deg=pu.get_dx10_map(surv='yr3',fwhm=np.radians(10.0))

hp.mollview m70f1qusm10deg[1],min=-2e-6,max=2e-6,title='DX10 full, RH1, 70 GHz Q, smth 10 degrees',unit='K'
savefig('qmap_full_rh1_70ghz_10deg.png')

hp.mollview m70f2qusm10deg[1],min=-2e-6,max=2e-6,title='DX10 full, RH2, 70 GHz Q, smth 10 degrees',unit='K'
savefig('qmap_full_rh2_70ghz_10deg.png')

hp.mollview m70f1qusm10deg[2],min=-2e-6,max=2e-6,title='DX10 full, RH1, 70 GHz U, smth 10 degrees',unit='K'
savefig('umap_full_rh1_70ghz_10deg.png')

hp.mollview m70f2qusm10deg[2],min=-2e-6,max=2e-6,title='DX10 full, RH2, 70 GHz U, smth 10 degrees',unit='K'
savefig('umap_full_rh2_70ghz_10deg.png')


m70y1ebsm1rad=pu.get_dx10_map(eb=True,surv='yr1',fwhm=1.0)
m70y2ebsm1rad=pu.get_dx10_map(eb=True,surv='yr2',fwhm=1.0)
m70y3ebsm1rad=pu.get_dx10_map(eb=True,surv='yr3',fwhm=1.0)

m70y1qusm1rad=pu.get_dx10_map(surv='yr1',fwhm=1.0)
m70y2qusm1rad=pu.get_dx10_map(surv='yr2',fwhm=1.0)
m70y3qusm1rad=pu.get_dx10_map(surv='yr3',fwhm=1.0)



m70f1ebsm=pu.get_dx10_map(eb=True,surv='full',half=1,fwhm=np.radians(10.0))
m70f2ebsm=pu.get_dx10_map(eb=True,surv='full',half=2,fwhm=np.radians(10.0))

m100f1ebsm=pu.get_dx10_map(chan='100',eb=True,surv='full',half=1,fwhm=np.radians(10.0))
m100f2ebsm=pu.get_dx10_map(chan='100',eb=True,surv='full',half=2,fwhm=np.radians(10.0))

m143f1ebsm=pu.get_dx10_map(chan='143',eb=True,surv='full',half=1,fwhm=np.radians(10.0))
m143f2ebsm=pu.get_dx10_map(chan='143',eb=True,surv='full',half=2,fwhm=np.radians(10.0))

m70f1qusm=pu.get_dx10_map(surv='full',half=1,fwhm=np.radians(10.0))
m70f2qusm=pu.get_dx10_map(surv='full',half=2,fwhm=np.radians(10.0))

m100f1qusm=pu.get_dx10_map(chan='100',surv='full',half=1,fwhm=np.radians(10.0))
m100f2qusm=pu.get_dx10_map(chan='100',surv='full',half=2,fwhm=np.radians(10.0))

m143f1qusm=pu.get_dx10_map(chan='143',surv='full',half=1,fwhm=np.radians(10.0))
m143f2qusm=pu.get_dx10_map(chan='143',surv='full',half=2,fwhm=np.radians(10.0))

hp.mollview m70f1qusm[1],min=-2e-7,max=2e-7,title='DX10 full, RH1, 70 GHz Q, smth 10 deg',unit='K'
savefig('qmap_full_rh1_70ghz_10deg.png')

hp.mollview m70f2qusm[1],min=-2e-7,max=2e-7,title='DX10 full, RH2, 70 GHz Q, smth 10 deg',unit='K'
savefig('qmap_full_rh2_70ghz_10deg.png')

hp.mollview m70f1qusm[2],min=-2e-7,max=2e-7,title='DX10 full, RH1, 70 GHz U, smth 10 deg',unit='K'
savefig('umap_full_rh1_70ghz_10deg.png')

hp.mollview m70f2qusm[2],min=-2e-7,max=2e-,title='DX10 full, RH2, 70 GHz U, smth 10 deg',unit='K'
savefig('umap_full_rh2_70ghz_10deg.png')


hp.mollview m100f1qusm[1],min=-2e-7,max=2e-7,title='DX10 full, RH1, 100 GHz Q, smth 10 deg',unit='K'
savefig('qmap_full_rh1_100ghz_10 deg.png')

hp.mollview m100f2qusm[1],min=-2e-7,max=2e-7,title='DX10 full, RH2, 100 GHz Q, smth 10 deg',unit='K'
savefig('qmap_full_rh2_100ghz_10deg.png')

hp.mollview m100f1qusm[2],min=-2e-7,max=2e-7,title='DX10 full, RH1, 100 GHz U, smth 10 deg',unit='K'
savefig('umap_full_rh1_100ghz_10deg.png')

hp.mollview m100f2qusm[2],min=-2e-7,max=2e-,title='DX10 full, RH2, 100 GHz U, smth 10 deg',unit='K'
savefig('umap_full_rh2_100ghz_10deg.png')


hp.mollview m143f1qusm1rad[1],min=-2e-7,max=2e-7,title='DX10 full, RH1, 143 GHz Q, smth 10 deg',unit='K'
savefig('qmap_full_rh1_143ghz_1rad.png')

hp.mollview m143f2qusm1rad[1],min=-2e-7,max=2e-7,title='DX10 full, RH2, 143 GHz Q, smth 10 deg',unit='K'
savefig('qmap_full_rh2_143ghz_1rad.png')

hp.mollview m143f1qusm1rad[2],min=-2e-7,max=2e-7,title='DX10 full, RH1, 143 GHz U, smth 10 deg',unit='K'
savefig('umap_full_rh1_143ghz_1rad.png')

hp.mollview m143f2qusm1rad[2],min=-2e-7,max=2e-,title='DX10 full, RH2, 143 GHz U, smth 10 deg',unit='K'
savefig('umap_full_rh2_143ghz_1rad.png')

#QU maps
chans=['070','100','143']
halfs=[1,2]
fwhms=[10.,40.]
scales=[2e-6,1e-6]
for chan in chans:
    for half in halfs:
        for fwhm,scale in zip(fwhms,scales):
            mqu=pu.get_dx10_map(chan=chan,half=half,fwhm=np.radians(fwhm),surv='full')
            hp.mollview(mqu[1],min=-1*scale,max=scale,title='DX10 full, RH'+np.str(half)+' ,  '+np.str(chan)+ 'GHz Q, smth' + np.str(fwhm)+' deg',unit='K')
            fname='qmap_full_rh'+np.str(half)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
            savefig(fname)
            hp.mollview(mqu[2],max=scale,min=-1.*scale,title='DX10 full, RH'+np.str(half)+' ,  '+np.str(chan)+ 'GHz U, smth' + np.str(fwhm)+' deg',unit='K')
            fname='umap_full_rh'+np.str(half)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
            savefig(fname)
            
#full qu,EB MAPS more channels
chans=['044','070','100','143','217']
survs=['full']
fwhms=[10.,40.]
scales=[2e-6,1e-6]
surv='full'
for chan in chans:
    for fwhm,scale in zip(fwhms,scales):
        meb=pu.get_dx10_map(eb=True,chan=chan,surv=surv,fwhm=np.radians(fwhm))
        hp.mollview(meb[1],min=-1*scale,max=scale,title='DX10 ,'+ np.str(surv)+' ,  '+np.str(chan)+ 'GHz E, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='emap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)
        hp.mollview(meb[2],max=scale,min=-1.*scale,title='DX10,'+np.str(surv)+' ,  '+np.str(chan)+ 'GHz B, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='bmap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)

for chan in chans:
    for fwhm,scale in zip(fwhms,scales):
        mqu=pu.get_dx10_map(chan=chan,surv=surv,fwhm=np.radians(fwhm))
        hp.mollview(mqu[1],min=-1*scale,max=scale,title='DX10 ,'+ np.str(surv)+' ,  '+np.str(chan)+ 'GHz Q, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='qmap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)
        hp.mollview(mqu[2],max=scale,min=-1.*scale,title='DX10,'+np.str(surv)+' ,  '+np.str(chan)+ 'GHz U, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='umap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)

survs=['full','nominal']
fwhm=2.
scale=4e-6
for chan in chans:
    for surv in survs:
        mqu=pu.get_dx10_map(chan=chan,surv=surv,fwhm=np.radians(fwhm))
        hp.mollview(mqu[1],min=-1*scale,max=scale,title='DX10 ,'+ np.str(surv)+' ,  '+np.str(chan)+ 'GHz Q, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='qmap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)
        hp.mollview(mqu[2],max=scale,min=-1.*scale,title='DX10,'+np.str(surv)+' ,  '+np.str(chan)+ 'GHz U, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='umap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)
        meb=pu.get_dx10_map(eb=True,chan=chan,surv=surv,fwhm=np.radians(fwhm))
        hp.mollview(meb[1],min=-1*scale,max=scale,title='DX10 ,'+ np.str(surv)+' ,  '+np.str(chan)+ 'GHz E, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='emap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)
        hp.mollview(meb[2],max=scale,min=-1.*scale,title='DX10,'+np.str(surv)+' ,  '+np.str(chan)+ 'GHz B, smth ' + np.str(fwhm)+' deg',unit='K')
        fname='bmap_'+np.str(surv)+'_'+np.str(chan)+'ghz'+np.str(fwhm)+'deg.png'
        savefig(fname)



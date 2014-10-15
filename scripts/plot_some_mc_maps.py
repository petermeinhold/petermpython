#first just get the mask by running a normal map

m1=pu.get_dx10_map(surv='full',half=1)

hp.mollview m1ebsm[2],min=-2e-6,max=2e-6,title='DX10 Year 1, 70 GHz B map, smoothed 10 degrees',unit='K'
savefig('bmap_yr1_70ghz_10deg.png')


f70rh1=pyfits.open('/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/lfi/LFI_SkyMap_070_1024_DX10_full_ringhalf_1.fits')
f70rh2=pyfits.open('/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/lfi/LFI_SkyMap_070_1024_DX10_full_ringhalf_2.fits')

f100rh1=pyfits.open('/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/hfi/official/HFI_SkyMap_100_2048_DX10_full_ringhalf_1.fits')
f100rh2=pyfits.open('/global/project/projectdirs/planck/data/mission/DPC_maps/dx10/hfi/official/HFI_SkyMap_100_2048_DX10_full_ringhalf_2.fits')

qqcov70_rh1=f70rh1[1].data['QQ_cov']
qqcov70_rh2=f70rh2[1].data['QQ_cov']
qqcov100_rh1=f100rh1[1].data['QQ_cov']
qqcov100_rh2=f100rh2[1].data['QQ_cov']

mc70_rh1=randn(len(qqcov70_rh1))*sqrt(qqcov70_rh1)
mc70_rh2=randn(len(qqcov70_rh2))*sqrt(qqcov70_rh2)

mc100_rh1=randn(len(qqcov100_rh1))*sqrt(qqcov100_rh1)
mc100_rh2=randn(len(qqcov100_rh2))*sqrt(qqcov100_rh2)

#note these are in Nest ordering, healpy assumes ring, so lets change, also degrade to 1024, change to masked array

mc70_rh1r=hp.ma(hp.ud_grade(mc70_rh1,1024,order_in='NEST',order_out='RING'))
mc70_rh2r=hp.ma(hp.ud_grade(mc70_rh2,1024,order_in='NEST',order_out='RING'))
mc100_rh1r=hp.ma(hp.ud_grade(mc100_rh1,1024,order_in='NEST',order_out='RING'))
mc100_rh2r=hp.ma(hp.ud_grade(mc100_rh2,1024,order_in='NEST',order_out='RING'))




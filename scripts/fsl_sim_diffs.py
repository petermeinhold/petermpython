fsld27m0=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=0)
fsld27m1=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=1)
fsld27m2=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=2)
fsld27m3=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=3)
fsld27m4=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=4)
fsld27m5=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=5)
fsld27m6=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='destriped',subband=6)

fsld27s0=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=0)
fsld27s1=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=1)
fsld27s2=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=2)
fsld27s3=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=3)
fsld27s4=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=4)
fsld27s5=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=5)
fsld27s6=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='destriped',subband=6)

fsld28m0=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=0)
fsld28m1=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=1)
fsld28m2=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=2)
fsld28m3=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=3)
fsld28m4=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=4)
fsld28m5=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=5)
fsld28m6=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='destriped',subband=6)

fsld28s0=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=0)
fsld28s1=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=1)
fsld28s2=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=2)
fsld28s3=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=3)
fsld28s4=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=4)
fsld28s5=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=5)
fsld28s6=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='destriped',subband=6)

hp.mollview fsld28m0,title='LFI28M f0 SS1-SS2',min=-20e-6,max=5e-6 
hp.mollview fsld28m1,title='LFI28M f1 SS1-SS2',min=-20e-6,max=5e-6 
hp.mollview fsld28m2,title='LFI28M f2 SS1-SS2',min=-20e-6,max=5e-6 
hp.mollview fsld28m3,title='LFI28M f3 SS1-SS2',min=-20e-6,max=5e-6 
hp.mollview fsld28m4,title='LFI28M f4 SS1-SS2',min=-20e-6,max=5e-6 
hp.mollview fsld28m5,title='LFI28M f5 SS1-SS2',min=-20e-6,max=5e-6 
hp.mollview fsld28m6,title='LFI28M f6 SS1-SS2',min=-20e-6,max=5e-6 



allsims=[]
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m0),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m1),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m2),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m3),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m4),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m5),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27m6),128,pess=True)))

allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s0),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s1),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s2),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s3),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s4),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s5),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld27s6),128,pess=True)))

allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m0),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m1),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m2),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m3),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m4),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m5),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28m6),128,pess=True)))

allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s0),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s1),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s2),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s3),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s4),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s5),128,pess=True)))
allsims.append(hp.ma(hp.ud_grade(maskmap(fsld28s6),128,pess=True)))





l28s12=pu.get_surv_diff_chan(chan1='28S',surv1=1,chan2='28S',surv2=2,nside=128)
l28s34=pu.get_surv_diff_chan(chan1='28S',surv1=3,chan2='28S',surv2=4,nside=128)

l27s12=pu.get_surv_diff_chan(chan1='27S',surv1=1,chan2='27S',surv2=2,nside=128)
l27s34=pu.get_surv_diff_chan(chan1='27S',surv1=3,chan2='27S',surv2=4,nside=128)

l28m12=pu.get_surv_diff_chan(chan1='28M',surv1=1,chan2='28M',surv2=2,nside=128)
l28m34=pu.get_surv_diff_chan(chan1='28M',surv1=3,chan2='28M',surv2=4,nside=128)

l27m12=pu.get_surv_diff_chan(chan1='27M',surv1=1,chan2='27M',surv2=2,nside=128)
l27m34=pu.get_surv_diff_chan(chan1='27M',surv1=3,chan2='27M',surv2=4,nside=128)





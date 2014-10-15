fsld27m0=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=0)
fsld27m1=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=1)
fsld27m2=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=2)
fsld27m3=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=3)
fsld27m4=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=4)
fsld27m5=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=5)
fsld27m6=pu.get_surv_diff_fsl(chan=27,rad='M',maptype='binned',subband=6)

fsld27s0=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=0)
fsld27s1=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=1)
fsld27s2=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=2)
fsld27s3=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=3)
fsld27s4=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=4)
fsld27s5=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=5)
fsld27s6=pu.get_surv_diff_fsl(chan=27,rad='S',maptype='binned',subband=6)

fsld28m0=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=0)
fsld28m1=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=1)
fsld28m2=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=2)
fsld28m3=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=3)
fsld28m4=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=4)
fsld28m5=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=5)
fsld28m6=pu.get_surv_diff_fsl(chan=28,rad='M',maptype='binned',subband=6)

fsld28s0=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=0)
fsld28s1=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=1)
fsld28s2=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=2)
fsld28s3=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=3)
fsld28s4=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=4)
fsld28s5=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=5)
fsld28s6=pu.get_surv_diff_fsl(chan=28,rad='S',maptype='binned',subband=6)


allsimsbin=[]
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m0),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m1),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m2),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m3),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m4),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m5),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27m6),128,pess=True)))

allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s0),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s1),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s2),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s3),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s4),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s5),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld27s6),128,pess=True)))

allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m0),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m1),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m2),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m3),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m4),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m5),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28m6),128,pess=True)))

allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s0),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s1),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s2),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s3),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s4),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s5),128,pess=True)))
allsimsbin.append(hp.ma(hp.ud_grade(maskmap(fsld28s6),128,pess=True)))







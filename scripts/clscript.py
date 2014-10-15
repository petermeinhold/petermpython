

dx930s1=pu.get_dx_map(surv='survey_1')
dx930s2=pu.get_dx_map(surv='survey_2')
dx930s3=pu.get_dx_map(surv='survey_3')
dx930s4=pu.get_dx_map(surv='survey_4')
dx930s5=pu.get_dx_map(surv='survey_5')

map=hp.ma(maskmap(bmap))
tmask=dx930s1.mask|dx930s2.mask|dx930s3.mask|dx930s4.mask|dx930s5.mask|bmap.mask

dx930s1.mask=tmask
dx930s2.mask=tmask
dx930s3.mask=tmask
dx930s4.mask=tmask
dx930s5.mask=tmask

cl12=hp.anafast(dx930s1.filled(),dx930s2.filled())
cl13=hp.anafast(dx930s1.filled(),dx930s3.filled())
cl14=hp.anafast(dx930s1.filled(),dx930s4.filled())
cl15=hp.anafast(dx930s1.filled(),dx930s5.filled())
cl23=hp.anafast(dx930s2.filled(),dx930s3.filled())
cl24=hp.anafast(dx930s2.filled(),dx930s4.filled())
cl25=hp.anafast(dx930s2.filled(),dx930s5.filled())
cl34=hp.anafast(dx930s3.filled(),dx930s4.filled())
cl35=hp.anafast(dx930s3.filled(),dx930s5.filled())
cl45=hp.anafast(dx930s4.filled(),dx930s5.filled())


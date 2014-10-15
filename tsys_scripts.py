import planck_util_prm as pu
import cPickle
import pyfits
from glob import glob
from scipy.stats import nanmean,nanstd,nanmedian
ucds30=pu.get_ucds(freq=30)
ucds44=pu.get_ucds(freq=44)
ucds70=pu.get_ucds(freq=70)
ucds=dict(ucds70,**ucds44)
ucds.update(ucds30)
ftswn=glob('C*tswn*7.fits')
fdcrs=glob('C*DCRS*7.fits')
ftswn.sort()
fdcrs.sort()
tswn30=pyfits.open(ftswn[0])
dcrs30=pyfits.open(fdcrs[0])
tswn44=pyfits.open(ftswn[1])
dcrs44=pyfits.open(fdcrs[1])
tswn70=pyfits.open(ftswn[2])
dcrs70=pyfits.open(fdcrs[2])
rcalist=[]
for r in range(18,29):
    for d in ['M','S']:
        rcalist.append('LFI'+str(r)+d)

for rca in rcalist[:13]:
    plt.figure(figsize=[12,5])   
    plt.plot (tswn70['pid'].data['pid'],tswn70[rca].data[rca],color='blue',label='Tsys gain'  )
    plt.plot (dcrs70['pid'].data['pid'],dcrs70[rca].data[rca],color='red',label='DCRS gain'  )
    plt.plot (ucds[rca]['pID'],ucds[rca]['nominal_gain'],color='black',label='UCDS nominal gain',lw=2)
    plt.xlabel('pid'),ylabel('Calibration, K/V'),title(rca)
    s=np.std(ucds[rca]['nominal_gain'])
    m=np.median(ucds[rca]['nominal_gain'])
    plt.ylim([m-5*s,m+5*s])
    leg=legend()
    grid()
    plt.show()
    

for rca in rcalist[13:19]:
    figure(figsize=[12,5])   
    plot (tswn44['pid'].data['pid'],tswn44[rca].data[rca],color='blue',label='Tsys gain'  )
    plot (dcrs44['pid'].data['pid'],dcrs44[rca].data[rca],color='red',label='DCRS gain'  )
    plot (ucds[rca]['pID'],ucds[rca]['nominal_gain'],color='black',label='UCDS nominal gain',lw=2)
    xlabel('pid'),ylabel('Calibration, K/V'),title(rca)
    leg=legend()
    

for rca in rcalist[19:]:
    figure(figsize=[12,5])   
    plot (tswn30['pid'].data['pid'],tswn30[rca].data[rca],color='blue',label='Tsys gain'  )
    plot (dcrs30['pid'].data['pid'],dcrs30[rca].data[rca],color='red',label='DCRS gain'  )
    plot (ucds[rca]['pID'],ucds[rca]['nominal_gain'],color='black',label='UCDS nominal gain',lw=2)
    xlabel('pid'),ylabel('Calibration, K/V'),title(rca)
    leg=legend()



odrange=arange(10)+945
chans=range(18,29)
alldata=[]
for rca in chans[:1]:
    ref=[]
    sky=[]
    obt=[]
    for od in odrange:
        o,s,r=pu.get_lfi_sky_ref(od,rca,'00')
        ref.append(r)
        sky.append(s)
        obt.append(o)
    obt=np.concatenate(obt)
    sky=np.concatenate(sky)
    ref=np.concatenate(ref)
    alldata.append({'obt':obt,'sky':sky,'ref':ref})
    
alldata=[]
for rca in chans:
    ref=[]
    sky=[]
    obt=[]
    for od in odrange:
        o,s,r=pu.get_lfi_sky_ref(od,rca,'01')
        ref.append(r)
        sky.append(s)
        obt.append(o)
    obt=np.concatenate(obt)
    sky=np.concatenate(sky)
    ref=np.concatenate(ref)
    alldata01.append({'obt':obt,'sky':sky,'ref':ref})
    
    
    
ref0=[]
sky0=[]
obt0=[]
for od in odrange:
    o,s,r=pu.get_lfi_sky_ref(od,rca,'00')
    ref0.append(r)
    sky0.append(s)
    obt0.append(o)
obt0=np.concatenate(obt0)
sky0=np.concatenate(sky0)
ref0=np.concatenate(ref0)


ref1=[]
sky1=[]
obt1=[]
for od in odrange:
    o,s,r=pu.get_lfi_sky_ref(od,rca,'01')
    ref1.append(r)
    sky1.append(s)
    obt1.append(o)
obt1=np.concatenate(obt1)
sky1=np.concatenate(sky1)
ref1=np.concatenate(ref1)


for key in sort(ucds70.keys()):
    thot=ucds70[key]['l14k']*.68
    thot[thot<1]=4.655*.68
    tcold=np.zeros(len(thot))+2.735*.51
    y0=(ucds70[key]['vref0'])/(ucds70[key]['vsky0']) 
    tsys0=(thot-y0*tcold)/(y0-1.)
    y1=(ucds70[key]['vref1'])/(ucds70[key]['vsky1']) 
    tsys1=(thot-y1*tcold)/(y1-1.)

    print key,' ref0  ',nanmedian(ucds70[key]['vref0']*ucds70[key]['gain']),nanmedian(tsys0)
    print key,' ref1  ',nanmedian(ucds70[key]['vref1']*ucds70[key]['gain']),nanmedian(tsys1)
    print key,' ref   ',nanmedian(ucds70[key]['vref']*ucds70[key]['gain'])
for key in sort(ucds44.keys()):
    thot=ucds44[key]['l14k']*.79
    thot[thot<1]=4.655*.79
    tcold=np.zeros(len(thot))+2.735*.66
    y0=(ucds44[key]['vref0'])/(ucds44[key]['vsky0']) 
    tsys0=(thot-y0*tcold)/(y0-1.)
    y1=(ucds44[key]['vref1'])/(ucds44[key]['vsky1']) 
    tsys1=(thot-y1*tcold)/(y1-1.)

    print key,' ref0  ',nanmedian(ucds44[key]['vref0']*ucds44[key]['gain']),nanmedian(tsys0)
    print key,' ref1  ',nanmedian(ucds44[key]['vref1']*ucds44[key]['gain']),nanmedian(tsys1)
    print key,' ref   ',nanmedian(ucds44[key]['vref']*ucds44[key]['gain'])
for key in sort(ucds30.keys()):
    thot=ucds30[key]['l14k']*.85
    thot[thot<1]=4.655*.85
    tcold=np.zeros(len(thot))+2.735*.76
    y0=(ucds30[key]['vref0'])/(ucds30[key]['vsky0']) 
    tsys0=(thot-y0*tcold)/(y0-1.)
    y1=(ucds30[key]['vref1'])/(ucds30[key]['vsky1']) 
    tsys1=(thot-y1*tcold)/(y1-1.)

    print key,' ref0  ',nanmedian(ucds30[key]['vref0']*ucds30[key]['gain']),nanmedian(tsys0)
    print key,' ref1  ',nanmedian(ucds30[key]['vref1']*ucds30[key]['gain']),nanmedian(tsys1)
    print key,' ref   ',nanmedian(ucds30[key]['vref']*ucds30[key]['gain'])

key='LFI22S'
offs=arange(11.)/10. -.5
for voff in offs:
    thot=ucds70[key]['l14k']*.68
    thot[thot<1]=4.655*.68
    tcold=np.zeros(len(thot))+2.735*.51
    y0=(ucds70[key]['vref0']-voff)/(ucds70[key]['vsky0']-voff) 
    tsys0=(thot-y0*tcold)/(y0-1.)
    y1=(ucds70[key]['vref1']-voff)/(ucds70[key]['vsky1']-voff) 
    tsys1=(thot-y1*tcold)/(y1-1.)
    print voff, key,' ref0  ',nanmedian((ucds70[key]['vref0']-voff)*ucds70[key]['gain']),nanmedian(tsys0)
    print voff, key,' ref1  ',nanmedian((ucds70[key]['vref1']-voff)*ucds70[key]['gain']),nanmedian(tsys1)
    print voff, key,' ref   ',nanmedian((ucds70[key]['vref']-voff)*ucds70[key]['gain'])    


rca=22
odrange=arange(10)+945
sref=[]
ssky=[]
sobt=[]
for od in odrange:
    o,s,r=pu.get_lfi_sky_ref(od,rca,'10')
    sref.append(r)
    ssky.append(s)
    sobt.append(o)
sobt=np.concatenate(sobt)
ssky=np.concatenate(ssky)
sref=np.concatenate(sref)
mref=[]
msky=[]
mobt=[]
for od in odrange:
    o,s,r=pu.get_lfi_sky_ref(od,rca,'11')
    mref.append(r)
    msky.append(s)
    mobt.append(o)
mobt=np.concatenate(mobt)
msky=np.concatenate(msky)
mref=np.concatenate(mref)




fracs=arange(11)/50.+.80
psds2=[]
for f in fracs:
    psds2.append(pu.nps((sref-rr*f*mref),78.5,minfreq=.1))
    

    
    
obt70=pu.rebin(alldata00[0]['obt'][:68000000],[1000],ratio=True)
obt44=pu.rebin(alldata00[6]['obt'][:40000000],[1000],ratio=True)
obt30=pu.rebin(alldata00[9]['obt'][:28000000],[1000],ratio=True)

for i,rca in enumerate(alldata00):
    if i < 6:
        r=pu.rebin(alldata00[i]['ref'][:68000000],[1000],ratio=True)
        plot (obt70,r-np.mean(r),label='rca '+str(i+18))
    if (i > 5) and ( i < 9):
        r=pu.rebin(alldata00[i]['ref'][:40000000],[1000],ratio=True)
        plot (obt44,r-np.mean(r),label='rca '+str(i+18))
    if i > 8:
        r=pu.rebin(alldata00[i]['ref'][:28000000],[1000],ratio=True)
        plot (obt30,r-np.mean(r),label='rca '+str(i+18))
        
        
        
#section to include diode weights to compare diode 0 and 1 with appropriate 'nominal gain'
Channel name	,	Weight
weight_dict={
"LFI18M00" :	0.567304963, 
"LFI18M01" :	0.432695037,
"LFI18S10" :	0.387168785,
"LFI18S11" :	0.612831215,
"LFI19M00" :	0.502457723,
"LFI19M01" :	0.497542277,
"LFI19S10" :	0.55143474,
"LFI19S11" :	0.44856526,
"LFI20M00" :	0.523020094,
"LFI20M01" :	0.476979906,
"LFI20S10" :	0.476730576,
"LFI20S11" :	0.523269424,
"LFI21M00" :	0.500324722,
"LFI21M01" :	0.499675278,
"LFI21S10" :	0.563712153,
"LFI21S11" :	0.436287847,
"LFI22M00" :	0.536283158,
"LFI22M01" :	0.463716842,
"LFI22S10" :	0.553913461,
"LFI22S11" :	0.446086539,
"LFI23M00" :	0.508036034,
"LFI23M01" :	0.491963966,
"LFI23S10" :	0.36160661,
"LFI23S11" :	0.63839339,
"LFI24M00" :	0.602269189,
"LFI24M01" :	0.397730811,
"LFI24S10" :	0.456037835,
"LFI24S11" :	0.543962165,
"LFI25M00" :	0.482050606,
"LFI25M01" :	0.517949394,
"LFI25S10" :	0.369618239,
"LFI25S11" :	0.630381761,
"LFI26M00" :	0.593126369,
"LFI26M01" :	0.406873631,
"LFI26S10" :	0.424268188,
"LFI26S11" :	0.575731812,
"LFI27M00" :	0.519877701,
"LFI27M01" :	0.480122299,
"LFI27S10" :	0.484831449,
"LFI27S11" :	0.515168551,
"LFI28M00" :	0.553227696,
"LFI28M01" :	0.446772304,
"LFI28S10" :	0.467677355,
"LFI28S11" :	0.532322645
}

rj2planck4k={
"LFI18M00" :	0.51, 
"LFI18M01" :	0.68,
"LFI18S10" :	0.68,
"LFI18S11" :	0.68,
"LFI19M00" :	0.68,
"LFI19M01" :	0.68,
"LFI19S10" :	0.68,
"LFI19S11" :	0.68,
"LFI20M00" :	0.68,
"LFI20M01" :	0.68,
"LFI20S10" :	0.68,
"LFI20S11" :	0.68,
"LFI21M00" :	0.68,
"LFI21M01" :	0.68,
"LFI21S10" :	0.68,
"LFI21S11" :	0.68,
"LFI22M00" :	0.68,
"LFI22M01" :	0.68,
"LFI22S10" :	0.68,
"LFI22S11" :	0.68,
"LFI23M00" :	0.68,
"LFI23M01" :	0.68,
"LFI23S10" :	0.68,
"LFI23S11" :	0.68,
"LFI24M00" :	0.79,
"LFI24M01" :	0.79,
"LFI24S10" :	0.79,
"LFI24S11" :	0.79,
"LFI25M00" :	0.79,
"LFI25M01" :	0.79,
"LFI25S10" :	0.79,
"LFI25S11" :	0.79,
"LFI26M00" :	0.79,
"LFI26M01" :	0.79,
"LFI26S10" :	0.79,
"LFI26S11" :	0.79,
"LFI27M00" :	0.85,
"LFI27M01" :	0.85,
"LFI27S10" :	0.85,
"LFI27S11" :	0.85,
"LFI28M00" :	0.85,
"LFI28M01" :	0.85,
"LFI28S10" :	0.85,
"LFI28S11" :	0.85
}

rj2plancksky={
"LFI18M00" :	0.51, 
"LFI18M01" :	0.51,
"LFI18S10" :	0.51,
"LFI18S11" :	0.51,
"LFI19M00" :	0.51,
"LFI19M01" :	0.51,
"LFI19S10" :	0.51,
"LFI19S11" :	0.51,
"LFI20M00" :	0.51,
"LFI20M01" :	0.51,
"LFI20S10" :	0.51,
"LFI20S11" :	0.51,
"LFI21M00" :	0.51,
"LFI21M01" :	0.51,
"LFI21S10" :	0.51,
"LFI21S11" :	0.51,
"LFI22M00" :	0.51,
"LFI22M01" :	0.51,
"LFI22S10" :	0.51,
"LFI22S11" :	0.51,
"LFI23M00" :	0.51,
"LFI23M01" :	0.51,
"LFI23S10" :	0.51,
"LFI23S11" :	0.51,
"LFI24M00" :	0.66,
"LFI24M01" :	0.66,
"LFI24S10" :	0.66,
"LFI24S11" :	0.66,
"LFI25M00" :	0.66,
"LFI25M01" :	0.66,
"LFI25S10" :	0.66,
"LFI25S11" :	0.66,
"LFI26M00" :	0.66,
"LFI26M01" :	0.66,
"LFI26S10" :	0.66,
"LFI26S11" :	0.66,
"LFI27M00" :	0.76,
"LFI27M01" :	0.76,
"LFI27S10" :	0.76,
"LFI27S11" :	0.76,
"LFI28M00" :	0.76,
"LFI28M01" :	0.76,
"LFI28S10" :	0.76,
"LFI28S11" :	0.76
}

#From villa, A&A prelaunch calibration paper
rca_tsys={"LFI18M00" :   36.0 ,	
"LFI18M01" :36.1,	
"LFI18S10" :33.9,	
"LFI18S11" :35.1,
"LFI19M00" :33.1,
"LFI19M01" :31.5,	
"LFI19S10" :32.2,	
"LFI19S11" :33.6,
"LFI20M00" :35.2,
"LFI20M01" :34.2,	
"LFI20S10" :36.9,	
"LFI20S11" :35.0,
"LFI21M00" :27.3,
"LFI21M01" :28.4,	
"LFI21S10" :34.4,	
"LFI21S11" :36.4,
"LFI22M00" :30.9,
"LFI22M01" :30.3,	
"LFI22S10" :30.3,	
"LFI22S11" :31.8,
"LFI23M00" :35.9,
"LFI23M01" :34.1,	
"LFI23S10" :33.9,	
"LFI23S11" :31.1,
"LFI24M00" :15.5,
"LFI24M01" :15.3,	
"LFI24S10" :15.8,	
"LFI24S11" :15.8,
"LFI25M00" :17.5,
"LFI25M01" :17.9,	
"LFI25S10" :18.6,	
"LFI25S11" :18.4,
"LFI26M00" :18.4,
"LFI26M01" :17.4,	
"LFI26S10" :16.8,	
"LFI26S11" :16.5,
"LFI27M00" :12.1,
"LFI27M01" :11.9,	
"LFI27S10" :13.0,	
"LFI27S11" :12.5,
"LFI28M00" :10.6,
"LFI28M01" :10.3,	
"LFI28S10" :9.9,
"LFI28S11" :9.8}


ucds30=pu.get_ucds(freq=30)
ucds44=pu.get_ucds(freq=44)
ucds70=pu.get_ucds(freq=70)
ucds=dict(ucds70,**ucds44)
ucds.update(ucds30)

#key='LFI22S'
#best choices so far offs=-0.0015
#                    offr=-0.03
offs=0.0
offr=0.0

#offsets=-.035+arange(21)/2000.
#print "Vrefoffset, Vsky offset, Tsys sky, Tsys ref, Tsys Y"
#for offr in offsets:
#imsky=np.zeros([100,100],dtype=np.float)
#imref=np.zeros([100,100],dtype=np.float)
#imy=np.zeros([100,100],dtype=np.float)
#for i in range(100):
#    offr=-.03+(i-50)/1000.
#    for j in range(100):
#        offs=-.0015+(j-50)/1000.

keys=ucds.keys()
keys.sort()
figure()
for key in keys:
    if key[-1]=='S':
        diode='10'
    if key[-1]=='M':
        diode='00'
    y=(ucds[key]['vref']-offr)/(ucds[key]['vsky']-offs)
    t4k=ucds[key]['l14k']
    t4k[t4k==0]=4.655
    thot=rj2planck4k[key+diode]*t4k
    tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(y)))
    tsys_ref=(ucds[key]['vref']-offr)*ucds[key]['nominal_gain']-thot
    tsys_sky=(ucds[key]['vsky']-offs)*ucds[key]['nominal_gain']-tcold
    tsys_y=(thot- y*tcold)/(y-1)
    imsky[i,j]=nanmedian(tsys_sky)
    imref[i,j]=nanmedian(tsys_ref)
    imy[i,j]=nanmedian(tsys_y)

    hold(False)
    plt.plot (ucds[key]['pID'],tsys_ref-nanmean(tsys_ref),label='Tsys_ref')
    plt.hold(True)

    plt.plot (ucds[key]['pID'],tsys_sky-nanmean(tsys_sky),label='Tsys_sky')
    plt.plot (ucds[key]['pID'],tsys_y-nanmean(tsys_y),label='Tsys_Y')
    plt.plot (ucds[key]['pID'],t4k-nanmean(t4k),label='l14k 4K temp')
    plt.plot (ucds[key]['pID'],ucds[key]['ts6r']-nanmean(ucds[key]['ts6r']),label='ts6R FP temp')
    plt.plot (ucds[key]['pID'],ucds[key]['lbem1']-nanmean(ucds[key]['lbem1']),label='lbem1 bem temp')
    plt.xlabel('pID'),plt.ylabel('Tsys (mean subtracted), K'),plt.title(key+' UCDS based Tsys estimates, Sky and ref use nominal_gain (DDX9)')
    plt.grid()
    plt.ylim(-3,3)
    plt.leg=legend()
    pu.thicklegendlines(leg)
    plt.savefig('tsys_no_offs'+key+'.png')
    
    
key='LFI22S'    
if key[-1]=='S':
    diode='10'
if key[-1]=='M':
    diode='00'
y=(ucds[key]['vref']-offr)/(ucds[key]['vsky']-offs)
t4k=ucds[key]['l14k']
t4k[t4k==0]=4.655
thot=rj2planck4k[key+diode]*t4k
tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(y)))
tsys_ref=(ucds[key]['vref']-offr)*ucds[key]['nominal_gain']-thot
tsys_sky=(ucds[key]['vsky']-offs)*ucds[key]['nominal_gain']-tcold
tsys_y=(thot- y*tcold)/(y-1)

hold(False)
plt.plot (ucds[key]['pID'],tsys_ref-nanmean(tsys_ref)+nanmean(tsys_y),label='Tsys_ref')
plt.hold(True)
plt.plot (ucds[key]['pID'],tsys_sky-nanmean(tsys_sky)+nanmean(tsys_y),label='Tsys_sky')
plt.plot (ucds[key]['pID'],tsys_y-nanmean(tsys_y)+nanmean(tsys_y),label='Tsys_Y')
plt.plot (ucds[key]['pID'],t4k-nanmean(t4k)+nanmean(tsys_y),label='l14k 4K temp')
plt.plot (ucds[key]['pID'],ucds[key]['ts6r']-nanmean(ucds[key]['ts6r'])+nanmean(tsys_y),label='ts6R FP temp')
plt.plot (ucds[key]['pID'],ucds[key]['lbem1']-nanmean(ucds[key]['lbem1'])+nanmean(tsys_y),label='lbem1 bem temp')
plt.xlabel('pID'),plt.ylabel('Tsys (mean subtracted), K'),plt.title(key+' UCDS based Tsys estimates, Sky and ref use nominal_gain (DDX9)')
plt.grid()
plt.ylim(-5+nanmean(tsys_y),5+nanmean(tsys_y))
leg=plt.legend()
pu.thicklegendlines(leg)


hold(False)
plt.plot (ucds[key]['pID'],tsys_ref-nanmean(tsys_ref),label='Tsys_ref')
plt.hold(True)
plt.plot (ucds[key]['pID'],tsys_sky-nanmean(tsys_sky),label='Tsys_sky')
plt.plot (ucds[key]['pID'],tsys_y-nanmean(tsys_y),label='Tsys_Y')
plt.plot (ucds[key]['pID'],t4k-nanmean(t4k),label='l14k 4K temp')
plt.plot (ucds[key]['pID'],ucds[key]['ts6r']-nanmean(ucds[key]['ts6r']),label='ts6R FP temp')
plt.plot (ucds[key]['pID'],ucds[key]['lbem1']-nanmean(ucds[key]['lbem1']),label='lbem1 bem temp')
plt.xlabel('pID'),plt.ylabel('Tsys (mean subtracted), K'),plt.title(key+' UCDS based Tsys estimates, Sky and ref use nominal_gain (DDX9)')
plt.grid()
plt.ylim(-3,3)
leg=plt.legend()
pu.thicklegendlines(leg)
#    print offr,offs,nanmedian(tsys_ref), nanmedian(tsys_sky),nanmedian(tsys_y)


#use preflight Tsys estimates to find offsets:
keys=ucds.keys()
tsys_y_all={}
for key in keys:
    if key[-1]=='S':
        diode='10'
        tsys_rca=(rca_tsys[key+'10'] + rca_tsys[key+'11'])/2.0
    if key[-1]=='M':
        diode='00'
        tsys_rca=(rca_tsys[key+'00'] + rca_tsys[key+'01'])/2.0
    t4k=ucds[key]['l14k']
    t4k[t4k==0]=4.655
    thot=rj2planck4k[key+diode]*t4k
    tcold=(2.735*rj2plancksky[key+diode]+np.zeros(len(ucds[key])))
    offr=nanmedian(ucds[key]['vref']-(tsys_rca+thot)/ucds[key]['nominal_gain'])
    offs=nanmedian(ucds[key]['vsky']-(tsys_rca+tcold)/ucds[key]['nominal_gain'])
    tsys_ref=(ucds[key]['vref']-offr)*ucds[key]['nominal_gain']-thot
    tsys_sky=(ucds[key]['vsky']-offs)*ucds[key]['nominal_gain']-tcold
    y=(ucds[key]['vref']-offr)/(ucds[key]['vsky']-offs)
    tsys_y=(thot- y*tcold)/(y-1)
    tsys_y_all.update({key:{"pID":ucds[key]["pID"],"tsys_y":tsys_y}})
    hold(False)
    plt.plot (ucds[key]['pID'],tsys_ref,label='Tsys_ref, offr= '+str(offr)[:6])
    plt.hold(True)
    plt.plot (ucds[key]['pID'],tsys_sky,label='Tsys_sky, offs='+str(offs)[:6])
    plt.plot (ucds[key]['pID'],tsys_y,label='Tsys_Y')`
    plt.plot (ucds[key]['pID'],t4k-nanmean(t4k)+tsys_rca,label='l14k 4K temp')
    plt.plot (ucds[key]['pID'],ucds[key]['ts6r']-nanmean(ucds[key]['ts6r'])+tsys_rca,label='ts6R FP temp')
    plt.plot (ucds[key]['pID'],ucds[key]['lbem1']-nanmean(ucds[key]['lbem1'])+tsys_rca,label='lbem1 bem temp')
    plt.xlabel('pID'),plt.ylabel('Tsys (mean subtracted), K'),plt.title(key+' UCDS based Tsys estimates, Sky and ref use nominal_gain (DDX9)')
    plt.grid()
    plt.ylim(tsys_rca-3,3+tsys_rca)
    leg=plt.legend()
    pu.thicklegendlines(leg)
    plt.savefig('tsys_test_offs_fit_rca'+key+'.png')
f=open('tsys_y_all.pkl','wb')
cPickle.dump(tsys_y_all,f)
f.close()
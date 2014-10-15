m1922nomt=hp.ma(m1922nom[0])
m1922nomq=hp.ma(m1922nom[1])
m1922nomu=hp.ma(m1922nom[2])
m1922nomt.mask=psmask|wmask
m1922nomq.mask=psmask|wmask
m1922nomu.mask=psmask|wmask
m1922=[m1922nomt.filled(),m1922nomq.filled(),m1922nomu.filled()]

m1823nomt=hp.ma(m1823nom[0])
m1823nomq=hp.ma(m1823nom[1])
m1823nomu=hp.ma(m1823nom[2])
m1823nomt.mask=psmask|wmask
m1823nomq.mask=psmask|wmask
m1823nomu.mask=psmask|wmask
m1823=[m1823nomt.filled(),m1823nomq.filled(),m1823nomu.filled()]

m2021nomt=hp.ma(m2021nom[0])
m2021nomq=hp.ma(m2021nom[1])
m2021nomu=hp.ma(m2021nom[2])
m2021nomt.mask=psmask|wmask
m2021nomq.mask=psmask|wmask
m2021nomu.mask=psmask|wmask
m2021=[m2021nomt.filled(),m2021nomq.filled(),m2021nomu.filled()]



clabtt,clabee,clabbb,clabte,clabtb,clabeb=hp.anafast(m1922,m1823)
clactt,clacee,clacbb,clacte,clactb,claceb=hp.anafast(m1922,m2021)
clbctt,clbcee,clbcbb,clbcte,clbctb,clbceb=hp.anafast(m1823,m2021)

plot l,l*(l+1)*clabtt,label='lfi19_22 x LFI 18_23'
title('Dx9 nominal survey TT wmap+pt src mask pseudospectrum')
xlabel('Multipole')
ylabel('l(l+1)Cl')
plot l,l*(l+1)*clactt,label='lfi19_22 x LFI 20_21'
plot l,l*(l+1)*clbctt,label='lfi18_23 x LFI 20_21'

figure()
plot l,l*(l+1)*clabee,label='lfi19_22 x LFI 18_23'
title('Dx9 nominal survey EE wmap+pt src mask pseudospectrum')
xlabel('Multipole')
ylabel('l(l+1)Cl')
plot l,l*(l+1)*clacee,label='lfi19_22 x LFI 20_21'
plot l,l*(l+1)*clbcee,label='lfi18_23 x LFI 20_21'

figure()
plot l,l*(l+1)*clabte,label='lfi19_22 x LFI 18_23'
title('Dx9 nominal survey TE wmap+pt src mask pseudospectrum')
xlabel('Multipole')
ylabel('l(l+1)Cl')
plot l,l*(l+1)*clacte,label='lfi19_22 x LFI 20_21'
plot l,l*(l+1)*clbcte,label='lfi18_23 x LFI 20_21'


plot l,l*(l+1)*clabee,label='lfi19_22 x LFI 18_23'
title('Dx9 nominal survey EE wmap+pt src mask pseudospectrum')
xlabel('Multipole')
ylabel('l(l+1)Cl')
plot l,l*(l+1)*clacee,label='lfi19_22 x LFI 20_21'
plot l,l*(l+1)*clbcee,label='lfi18_23 x LFI 20_21'

claatt,claaee,claabb,claate,claatb,claaeb=hp.anafast(m1922)
clbbtt,clbbee,clbbbb,clbbte,clbbtb,clbbeb=hp.anafast(m1823)
clcctt,clccee,clccbb,clccte,clcctb,clcceb=hp.anafast(m2021)


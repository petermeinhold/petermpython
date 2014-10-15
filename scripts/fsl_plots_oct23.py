
from petermpython import planck_util as pu
import healpy as hp
from glob import glob

basefsl='/global/homes/p/peterm/gpeterm/sidelobes/ddx9/'
#fl70fsl=glob(basefsl+'70ghz/*full_1024.fits')
#fsl70tqu=hp.ma(hp.read_map(fl70fsl[0],[0,1,2]))

#cl70=hp.anafast(fsl70tqu)
l=arange(len(cl70[0]))
figure(figsize=(10,4))

plot l,1e12*cl70[0]*l*(l+1)/(2*np.pi),'b-',label='TT'
plot l,1e12*cl70[1]*l*(l+1)/(2*np.pi),'g-',label='EE'
plot l,1e12*cl70[2]*l*(l+1)/(2*np.pi),'r-',label='BB'

xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2') 
leg=legend()
pu.thicklegendlines(leg)
title('70 GHz monochromatic FSL model, Delta DX9 as input sky')
xlim[0,50]
ylim[0,.01]
grid()
savefig('plots/fsl_sims_cl_oct18/cl_fsl70.png')

basefsl='/global/homes/p/peterm/gpeterm/sidelobes/ddx9/'
fl30fsl=glob(basefsl+'30ghz/*full_1024.fits')
fsl30tqu=hp.ma(hp.read_map(fl30fsl[0],[0,1,2]))

#cl30=hp.anafast(fsl30tqu)
l=arange(len(cl30[0]))
figure(figsize=(10,4))
plot l,1e12*cl30[0]*l*(l+1)/(2*np.pi),'b-',label='TT'
plot l,1e12*cl30[1]*l*(l+1)/(2*np.pi),'g-',label='EE'
plot l,1e12*cl30[2]*l*(l+1)/(2*np.pi),'r-',label='BB'
xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2') 
leg=legend()
pu.thicklegendlines(leg)
title('30 GHz bandweighted FSL model, Delta DX9 as input sky')
xlim[0,50]
ylim[0,.25]
grid()
savefig('plots/fsl_sims_cl_oct18/cl_fsl30.png')

basefsl='/global/homes/p/peterm/gpeterm/sidelobes/ddx9/'
#fl44fsl=glob(basefsl+'44ghz/*full_1024.fits')
#fsl44tqu=hp.ma(hp.read_map(fl44fsl[0],[0,1,2]))

#cl44=hp.anafast(fsl44tqu)
l=arange(len(cl44[0]))
figure(figsize=(10,4))

plot l,1e12*cl44[0]*l*(l+1)/(2*np.pi),'b-',label='TT'
plot l,1e12*cl44[1]*l*(l+1)/(2*np.pi),'g-',label='EE'
plot l,1e12*cl44[2]*l*(l+1)/(2*np.pi),'r-',label='BB'

xlabel('l'),ylabel('l*(l+1)*Cl/(2*pi), microK^2') 
leg=legend()
pu.thicklegendlines(leg)

title('44 GHz monochromatic FSL model, Delta DX9 as input sky')
xlim[0,50]
ylim[0,.01]
grid()
savefig('plots/fsl_sims_cl_oct18/cl_fsl44.png')

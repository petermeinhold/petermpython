import numpy as np
import matplotlib.pyplot as plt


def getnoisebandwidth(ucds,channel,calnoise,bparam):
	"""Estimate Beta from noise and UCDS DC level and smoothed gain
	Assume calibrated noise is in RJ microK*sqrt(sec) Jan 2012 add nonlinear parameter correction"""
	k=(2.)**.5
	sci='SCI'+channel+'W'
	sg='SG'+channel
	sci_pid_sort=ucds[sci].data['PID'].ravel().argsort()
	sg_pid_sort=ucds[sg].data['PID'].ravel().argsort()
	scipid=ucds[sci].data['PID'][sci_pid_sort]
	sgpid=ucds[sg].data['PID'][sg_pid_sort]
	scisky=ucds[sci].data['MEANSKY'][sci_pid_sort]
	sggain=ucds[sg].data['GAIN'][sg_pid_sort]
	pidlist = np.arange(np.min(scipid)+150, np.max(scipid)-150,100)
	beta=np.zeros(len(pidlist),dtype=float)
        j=0
        for pidcen in pidlist :
		meangain=np.mean(sggain[abs(sgpid-pidcen) <50])
               	meandcsky=np.mean(scisky[abs(scipid-pidcen)<50])
		if ((meangain != 0) and (meandcsky !=0)):
		        beta[j]=((meangain*k*meandcsky/calnoise)/(1+bparam*meandcsky))**2.
	        j=j+1
	return(pidlist,beta)

def get_all_betas(ucds,calnoise,bparamdic):
       """run thru all chans, calnoise hardwired to published early flight results"""
       calnoise=np.array([453,412,515,492,523,552,403,499,435,473,449,480,440,380,394,374,456,399,272,297,306,280])*1.e-6

       beta=np.zeros((22,2),dtype=float)
       i=0
       for rad in np.arange(18,29):
	       horn='LFI'+str(rad)
               channelm=str(rad)+'0'
               channels=str(rad)+'1'
               print channels
               bparam_m=np.mean(bparamdic[horn][0:1])
               bparam_s=np.mean(bparamdic[horn][2:3])
	       print bparam_m,bparam_s
               pid,betam=getnoisebandwidth(ucds,channelm,calnoise[i*2],bparam_m)
               pid,betas=getnoisebandwidth(ucds,channels,calnoise[i*2+1],bparam_s)
  	       beta[2*i,0]=np.mean(betam[betam > 0])
               beta[2*i+1,0]=np.mean(betas[betas > 0])
               beta[2*i,1]=np.std(betam[betam > 0])
               beta[2*i+1,1]=np.std(betas[betas > 0])
               i=i+1
       return(beta)


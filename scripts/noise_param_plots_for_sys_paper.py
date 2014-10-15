from setup_matplotlib import *
from matplotlib.ticker import MaxNLocator
import cPickle
import numpy as np

fk_18m=np.genfromtxt('fk_marge/LFI18M.txt')
fk_24m=np.genfromtxt('fk_marge/LFI24M.txt')
fk_27m=np.genfromtxt('fk_marge/LFI27M.txt')

alpha_18m=np.genfromtxt('alpha_marge/LFI18M.txt')
alpha_24m=np.genfromtxt('alpha_marge/LFI24M.txt')
alpha_27m=np.genfromtxt('alpha_marge/LFI27M.txt')

wnlevel_18m=np.genfromtxt('wnlevel_marge/LFI18M.txt')
wnlevel_24m=np.genfromtxt('wnlevel_marge/LFI24M.txt')
wnlevel_27m=np.genfromtxt('wnlevel_marge/LFI27M.txt')
#for width in [18., 12., 8.8]:
for width in [8.8]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    
    plt.plot(wnlevel_18m[:,0]+2.5,wnlevel_18m[:,1],'+',label='LFI18M')
    plt.plot(wnlevel_24m[:,0]+2.5,wnlevel_24m[:,1],'o',label='LFI24M')
    plt.plot(wnlevel_27m[:,0]+2.5,wnlevel_27m[:,1],'x',label='LFI27M')
    plt.legend(frameon=False,loc=0)
    plt.grid()
    plt.xlim([80,600])
    # labels
    plt.xlabel(r"Days after launch"); plt.ylabel(r"White noise $\left[ \mu K_{cmb} \sqrt{Sec} \right]$")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # reduce white space around figure
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")
    plt.show()
    plt.savefig("wnlevel_%dmm.pdf" % int(width*10), bbox_inches='tight',pad_inches=0)
    plt.savefig("wnlevel_%dmm.eps" % int(width*10), bbox_inches='tight',pad_inches=0)
    
for width in [8.8]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    
    plt.plot(alpha_18m[:,0]+2.5,alpha_18m[:,1],'+',label='LFI18M')
    plt.plot(alpha_24m[:,0]+2.5,alpha_24m[:,1],'o',label='LFI24M')
    plt.plot(alpha_27m[:,0]+2.5,alpha_27m[:,1],'x',label='LFI27M')
    plt.legend(frameon=False,loc=0)
    plt.grid()
    plt.ylim([-1.5,-.85])
    plt.xlim([80,600])
    # labels
    plt.xlabel(r"Days after launch"); plt.ylabel(r"$\alpha$")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # reduce white space around figure
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")
    plt.show()
    plt.savefig("alpha_%dmm.pdf" % int(width*10), bbox_inches='tight',pad_inches=0)
    plt.savefig("alpha_%dmm.eps" % int(width*10), bbox_inches='tight',pad_inches=0)

for width in [8.8]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
    # this should be changed for making a panel of multiple figures
    ax = fig.add_subplot(111)
    
    plt.plot(fk_18m[:,0]+2.5,fk_18m[:,1],'+',label='LFI18M')
    plt.plot(fk_24m[:,0]+2.5,fk_24m[:,1],'o',label='LFI24M')
    plt.plot(fk_27m[:,0]+2.5,fk_27m[:,1],'x',label='LFI27M')
    plt.legend(frameon=False,loc=0)
    plt.grid()
    plt.ylim([0,200])
    plt.xlim([80,600])
    # labels
    plt.xlabel(r"Days after launch"); plt.ylabel(r"Knee Frequency, mHz")
    ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels
    # reduce ticks for small figures
    if width < 10:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    # reduce white space around figure
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    # set vertical y axis ticklables
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_rotation("vertical")
    plt.show()
    plt.savefig("fk_%dmm.pdf" % int(width*10), bbox_inches='tight',pad_inches=0)
    plt.savefig("fk_%dmm.eps" % int(width*10), bbox_inches='tight',pad_inches=0)
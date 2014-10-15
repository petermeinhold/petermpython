def get_rel_scales(psmap):
    """
    function to find linear fit params among adjacent
    phase binned rings in a pseudomap
    returns array of slopes and offsets
    row_i ~ fitpar1*row_(i+1) + fitpar2
    """
    shape(psmap)
    s=psmap.shape()
    nrows=s[0]
    fitpars1=[]
    fitpars2=[]
    for i in range(nrows-1):
        fp1,fp2=polyfit(psmap[i,:],psmap[i+1,:],1)
        fitpars1.append(fp1)
        fitpars2.append(fp2)
    fitpars1=np.array(fitpars1)
    fitpars2=np.array(fitpars2)
    return fitpars1,fitpars2

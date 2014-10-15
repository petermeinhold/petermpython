for chan in [30,44,70]:
    for surv1 in [1,2,3,4,5]:
        for surv2 in [1,2,3,4,5]:
            if surv2 > surv1:
                diff=pu.get_surv_diff(chan1=chan,chan2=chan,surv1=surv1,surv2=surv2,cal1='dx9',cal2='dx9')
                hp.mollview(diff,min=-25e-6,max=25e-6,title='DX9 '+ np.str(chan)+' GHz '+'Survey ' + np.str(surv1) + '- Survey '+ np.str(surv2))
                plt.savefig('dx9_'+np.str(chan)+'_s'+np.str(surv1)+'_s'+np.str(surv2)+'.png')


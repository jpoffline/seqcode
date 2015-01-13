# function to return correct labels

def mcmc_labels(ins):

    N = len(ins)

    dtoplot = []
    lab = []
    dloc = []
    fac = []
    for n in xrange(0,N):
        dtoplot.append(ins[n])
        fac.append(1.0)        

    for i in xrange(0,N):

    	lab.append('label')
    	ll = dtoplot[i];
    	if ll == 'phidot0':
    		lab[i] = '$\dot{\phi}_0$'
    		dloc.append(0)
    	if ll == 'm3':
    		lab[i] = '$m^3$'
    		dloc.append(1)
    	if ll == 'like':
    		lab[i] = 'L'
    		dloc.append(2)
    	if ll == 'Omk': 
    		lab[i] = '$10^{-3}\Omega_k h^2$'
    		dloc.append(3)
    		fac[i] = 1000.0
    	if ll == 'amax': 
    		lab[i] = '$a_{max}$'
    		dloc.append(4)
    	if ll == 'tmax': 
    		lab[i] = '$t_{max}$'
    		dloc.append(5)
    	if ll == 'tmaxfrac': 
    		lab[i] = '$t_{max} / t_0$'
    		dloc.append(6)
    	if ll == 'HIR': 
    		lab[i] = '$< R>$'
    		dloc.append(7)
            
    return (dloc, lab, fac)    
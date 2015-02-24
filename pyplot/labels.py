# function to return correct labels

def findloc(in_key, in_list):
    n = len(in_list)
    for i in xrange(0, n):
        if in_list[i] == in_key:
            return i

def read_priors(priorsfilename):
    # Open up priors file, containing labels of quanties actually MCMCd over
    priors_file = open(priorsfilename,'r');
    lines = priors_file.readlines();
    priors_file.close();

    legit_labels = []
    for line in lines:
    	if line.startswith('#') or not line.strip():
    		continue
    	else:
            p = line.split()
            to_put = p[1]
            if to_put == 'mass':
                to_put = 'm3'
            legit_labels.append(to_put)
            
    return legit_labels

def mcmc_labels(ins, priorsfilename):

    
    legit_labels = read_priors(priorsfilename)
    n_legit_labels = len(legit_labels)   
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
            dloc.append( findloc(ll, legit_labels) )
    	if ll == 'm3':
            lab[i] = '$m^3$'
            dloc.append( findloc(ll, legit_labels) ) 
        if ll == 'Omegamh2':
            lab[i] = '$\Omega_{\\rm m} h^2$'     
            dloc.append( findloc(ll, legit_labels) ) 
        if ll == 'Omegabh2':
            lab[i] = '$\Omega_{\\rm b} h^2$'     
            dloc.append( findloc(ll, legit_labels) ) 
    	if ll == 'Omegakh2': 
            lab[i] = '$10^{-3}\Omega_{\\rm k} h^2$'
            dloc.append( findloc(ll, legit_labels) ) 
            fac[i] = 1000.0
            

    	if ll == 'like':
    		lab[i] = 'L'
    		dloc.append(n_legit_labels)
    	if ll == 'Omk': 
    		lab[i] = '$10^{-3}\Omega_{\\rm k} h^2$'
    		dloc.append( n_legit_labels + 1 )
    		fac[i] = 1000.0
    	if ll == 'amax': 
    		lab[i] = '$a_{\\rm max} / a_0$'
    		dloc.append( n_legit_labels + 2 )
    	if ll == 'tmax': 
    		lab[i] = '$t_{\\rm max}$'
    		dloc.append( n_legit_labels + 3 )
    	if ll == 'tmaxfrac': 
    		lab[i] = '$t_{\\rm max} / t_0$'
    		dloc.append( n_legit_labels + 4 )
    	if ll == 'HIR': 
    		lab[i] = '$< R>$'
    		dloc.append( n_legit_labels + 5 )
                
    return (dloc, lab, fac)    
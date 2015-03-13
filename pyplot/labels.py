# function to return correct labels
import backend
import plottools

def GetFileNames(ID, chainID):
    data_DIR = '../' + backend.chainsDIR + '/' + ID + '/'
    priors_filename = data_DIR + backend.priors_file_name
    data_filename = data_DIR + chainID + '.dat'
    PLOTSDIR = backend.plotsDIR + '/' + ID + '/';
    plottools.check_dir_exists(PLOTSDIR)
    return (priors_filename, data_filename, PLOTSDIR)



def GetPriorsFileName(ID):
    return '../' + backend.chainsDIR + '/' + ID + '/' + backend.priors_file_name 

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
    
def read_priors_vals_name(priorsfilename):
    # Open up priors file, containing labels of quanties actually MCMCd over
    priors_file = open(priorsfilename,'r');
    lines = priors_file.readlines();
    priors_file.close();

    params_priors = ''
    for line in lines:
    	if line.startswith('#') or not line.strip():
    		continue
    	else:
            p = line.split()
            param_name = p[1]
            param_prior_lower = str(p[2])
            param_prior_upper = str(p[3])
            param_prior_sigma = str(p[4])

            params_priors += param_name + ' ' + param_prior_lower + ' ' + param_prior_upper + ' ' + param_prior_sigma + '\n'
            
    return params_priors    

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
            
    	if ll == 'phidot0':
            lab[i] = '$\dot{\phi}_0$'
            dloc.append( findloc(ll, legit_labels) )
            
    	if ll == 'm3':
            lab[i] = '$m^3$'
            dloc.append( findloc(ll, legit_labels) ) 
            

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
    
def getparams(priors_file_name):  

    all = read_priors(priors_file_name)
    all.append('Omk')
    all.append('amax')
    all.append('tmax')
    all.append('tmaxfrac')
    
    return all

    
    
def getallpairs(priors_file_name):
    (dloc, lab, fac) = mcmc_labels( getparams(priors_file_name) , priors_file_name)  
    return (dloc, lab)
    
def getLoc(param_name, priors_file_name):
    
    (dloc, lab, fac) = mcmc_labels( getparams(priors_file_name) , priors_file_name)
    all = getparams(priors_file_name)
    ID = False
    for i in xrange(0, len(lab)):
        if all[i] == param_name:
            ID = int(dloc[i])
            break
            
    return ID

def CheckDerived(toplot, priors_filename):
    isDERIVED = ''
    priors_MCMC_params = read_priors(priors_filename)
    for t in toplot:
        if t not in priors_MCMC_params:
            isDERIVED = backend.derived_param_prefix
            
    return isDERIVED        
    
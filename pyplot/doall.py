########################################################################
#
#   Driver to plot 1D & 2D marginalised distributions
#    obtained from MCMC
#   J Pearson, Nottingham, Jan 2015
#
#   Useage: from a bash terminal, run
#
#       > python doall.py RUNID
#
#   ...where RUNID is the name of the folder containing the chain
#      plots will be output into the directory 'plots/RUNID/'
#
########################################################################

import os, sys
import labels as lab
import hist2d_m
import hist1d
import backend

RUN_ID = sys.argv[1]
chainID = 'chain0001'

exceptions = ('tmax')

def main(RUN_ID, chainID):

    # Get the file name of the priors file
    priors_file_name = lab.GetPriorsFileName( RUN_ID )
    
    # Parameters to sweep over (to make plots of)
    all = lab.getparams(priors_file_name)

    if backend.verbose > -1:
        print ''
        print '<><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
        print 'Read in and compute 1D & 2D marginalised likelihoods'
        print '<><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
        print 'DATA: ', '../' + backend.chainsDIR + '/' + RUN_ID + '/'
        print 'PLOTS:', backend.plotsDIR + '/' + RUN_ID + '/'
        print ''
        print 'Params, priors, and sigmas used in MCMC'
        print lab.read_priors_vals_name(priors_file_name)
        params = 'All parameters, including derived ones:\n'  
        for param in all:
            params += param + ' '
        print params
        params = 'Ignored parameters: '
        for param in backend.exceptions:
            params += param + ' '
        print params    
        print ''
        if backend.mod_tmaxfrac:
            print 'tmax plot truncated between', backend.my_min_tmaxfrac, 'and', backend.my_max_tmaxfrac
            
    
    # Go get the relevant file names       
    filenames = lab.GetFileNames(RUN_ID, chainID)       
           
    # Run over all parameters to generate 1D & 2D marginalised likelihoods
    for i in xrange( 0 , len( all ) ):
        if all[i] not in backend.exceptions:
            if backend.verbose > 0:
                print '1D', all[i]
            hist1d.generate_histogram_1d(RUN_ID, [all[i]], chainID, filenames)
            for j in xrange( i + 1 , len( all ) ):
                if all[j] not in backend.exceptions:
                    if backend.verbose > 0:
                        print '2D', all[i], all[j]  
                    hist2d_m.generate_histogram_2d( RUN_ID, ( all[i], all[j] ) , chainID, filenames)


    print backend.done_message
    

main(RUN_ID, chainID)  
    
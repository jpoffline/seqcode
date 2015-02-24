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
#      plots will be output into the directory 'plots/'
#
########################################################################

import os, sys
import labels as lab
import hist2d_m
import hist1d

RUN_ID = sys.argv[1]

priors_file_name = '../chains/' + RUN_ID + '/priors.txt'

chainID = 'chain0001'


# Parameters to sweep over (to make plots of)
all = lab.read_priors(priors_file_name)
all.append('Omk')
all.append('amax')
all.append('tmaxfrac')

string = 'Parameters: '
for param in all:
    string += param + '\t'
print string

PY = 'python2.7'
print 'Read in and compute likelihoods'

for i in xrange(0, len(all)):
    print '1D', all[i]
    hist1d.generate_histogram_1d(RUN_ID, [all[i]], chainID)
    for j in xrange(i+1,len(all)):
        print '2D', all[i], all[j]  
        hist2d_m.generate_histogram_2d( RUN_ID, ( all[i], all[j] ) , chainID)

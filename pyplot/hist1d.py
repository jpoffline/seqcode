import plottools as pts
import labels
import pylab as P
import sys
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

# How many bins to use in the histogram?
nbins = 20

# Font-size for axes labels
fig_fontsize = 13

# Get "what to plot" from command line inputs
toplot = []
ID = 'run'
if len(sys.argv) > 1:
    ID = sys.argv[1]
    toplot.append(sys.argv[2])
else:
    # If nothing came through, use a default
    toplot.append('Omk')

data_DIR = '../chains/' + ID + '/'
priors_filename = data_DIR + 'priors.txt'
data_filename = data_DIR + 'chain0001.dat'
PLOTSDIR = 'plots/';

# Construct file name for outputted figure
out_fig_filename = PLOTSDIR + ID + '_' + toplot[0] +'_1dlike.pdf'


# Sort out data location & labels
dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

# Get the data columns
data = pts.data_1d(data_filename, dloc[0], fac[0])

# Obtain histogram
pts.plot_1dhist(data, nbins,axes_labels,out_fig_filename, fig_fontsize)
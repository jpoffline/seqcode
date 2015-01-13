import plottools as pts
import labels
import pylab as P
import sys
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
data_filename = '../chains/run2/chain0001.dat'
PLOTSDIR = 'plots/';
# How many bins to use in the histogram?
nbins = 20

# Font-size for axes labels
fig_fontsize = 13

# Get "what to plot" from command line inputs
toplot = []
if len(sys.argv) > 1:
    toplot.append(sys.argv[1])
else:
    # If nothing came through, use a default
    toplot.append('Omk')

# Construct file name for outputted figure
out_fig_filename = PLOTSDIR+toplot[0]+'_1dlike.pdf'


# Sort out data location & labels
dloc, axes_labels, fac = labels.mcmc_labels(toplot)

# Get the data columns
data = pts.data_1d(data_filename, dloc[0], fac[0])

# Obtain histogram
pts.plot_1dhist(data, nbins,axes_labels,out_fig_filename, fig_fontsize)
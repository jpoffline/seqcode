
import plottools as pts
import labels

from pylab import *

# Chain file name
data_filename = '../chains/run2/chain0001.dat'
PLOTSDIR = 'plots/';
# How many bins for the histogram
nbins = 40

# Font-size for the axes labels
fig_fontsize = 13

# Get "what to plot" from the command line
toplot = []
if len(sys.argv) > 1:
    toplot.append(sys.argv[1])
    toplot.append(sys.argv[2])
else:
    # Keep a set of defaults
    toplot.append('Omk')
    toplot.append('m3')

# Construct a file name
out_fig_filename = PLOTSDIR + toplot[0] + '_' + toplot[1] + '_2dlike.pdf'

# Sort out data location & labels
dloc, axes_labels, fac = labels.mcmc_labels(toplot)

# Get the data columns
data = pts.data_column(data_filename, dloc, fac)

# Obtain histogram
pts.plot_2dhist(data, axes_labels, nbins, fig_fontsize, out_fig_filename)
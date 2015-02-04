
import plottools as pts
import labels
from pylab import *


def generate_histogram_2d(ID, toplot):
    # How many bins for the histogram
    nbins = 40

    # Font-size for the axes labels
    fig_fontsize = 13

    # Get "what to plot" from the command line
  #  toplot = []
  #  if len(sys.argv) > 1:
  #      ID = sys.argv[1]
  #      toplot.append(sys.argv[2])
  #      toplot.append(sys.argv[3])
  #  else:
  #      # Keep a set of defaults
  #      toplot.append('Omk')
  #      toplot.append('m3')

    # Chain file name
    data_DIR = '../chains/'+ ID + '/'
    priors_filename = data_DIR + 'priors.txt'
    data_filename = data_DIR + 'chain0001.dat'
    PLOTSDIR = 'plots/';

    # Construct a file name
    out_fig_filename = PLOTSDIR + ID + '_' + toplot[0] + '_' + toplot[1] + '_2dlike.pdf'

    # Sort out data location & labels
    dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

    # Get the data columns
    data = pts.data_column(data_filename, dloc, fac)

    # Obtain histogram
    pts.plot_2dhist(data, axes_labels, nbins, fig_fontsize, out_fig_filename)
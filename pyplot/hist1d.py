import plottools as pts
import labels
import pylab as P
import sys
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

import backend

def generate_histogram_1d(ID, toplot, chainID):

    print toplot

    # How many bins to use in the histogram?
    nbins = backend.nbins_1d

    #chainID = 'chain0001'

    # Font-size for axes labels
    fig_fontsize = backend.fig_fontsize

    data_DIR = '../chains/' + ID + '/'
    priors_filename = data_DIR + backend.priors_file_name
    data_filename = data_DIR + chainID + '.dat'
    PLOTSDIR = 'plots/' + ID + '/';
    pts.check_dir_exists(PLOTSDIR)

    # Construct file name for outputted figure
    out_fig_filename = PLOTSDIR + ID + '_' + toplot[0] + backend.like_1d_fn

    print len(toplot)
    print out_fig_filename

    # Sort out data location & labels
    dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

    # Get the data columns
    data = pts.data_1d(data_filename, dloc[0], fac[0])

    # Obtain histogram
    pts.plot_1dhist(data, nbins,axes_labels,out_fig_filename, fig_fontsize)

import plottools as pts
import labels
from pylab import *
import backend

def generate_histogram_2d(ID, toplot, chainID):
    # How many bins for the histogram
    nbins = backend.nbins_2d

    # Font-size for the axes labels
    fig_fontsize = backend.fig_fontsize

    # Chain file name
    data_DIR = '../chains/'+ ID + '/'
    priors_filename = data_DIR + backend.priors_file_name
    data_filename = data_DIR + chainID + '.dat'
    PLOTSDIR = 'plots/' + ID + '/';
    pts.check_dir_exists(PLOTSDIR)

    # Construct a file name
    out_fig_filename = PLOTSDIR + ID + '_' + toplot[0] + '_' + toplot[1] + backend.like_2d_fn

    # Sort out data location & labels
    dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

    # Get the data columns
    data = pts.data_column(data_filename, dloc, fac)

    # Obtain histogram
    pts.plot_2dhist(data, axes_labels, nbins, fig_fontsize, out_fig_filename)
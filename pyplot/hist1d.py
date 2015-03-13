import plottools as pts
import labels
import pylab as P
import sys
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
import backend

def generate_histogram_1d(ID, toplot, chainID, filenames):

    # Unpack filenames from input
    (priors_filename, data_filename, PLOTSDIR) = filenames
    
    # Check whether any of the given parameters are "derived",
    # since all plots containing derived parameters are appended with a character.
    isDERIVED = labels.CheckDerived(toplot, priors_filename)

    # Construct file name for outputted figure
    out_fig_filename = PLOTSDIR + isDERIVED + ID + '_' + toplot[0] + backend.like_1d_fn

    # How many bins to use in the histogram?
    nbins = backend.nbins_1d
    
    # Font-size for axes labels
    fig_fontsize = backend.fig_fontsize    

    if backend.verbose > 1:
        print 'OUTPUT TO:', out_fig_filename

    # Sort out data location & labels
    dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

    # Get the data columns
    data = pts.data_1d(data_filename, dloc[0], fac[0])

    
    for t in toplot:
        if t == 'tmaxfrac':
            if backend.mod_tmaxfrac:
                data = pts.mod_data_tmax(data, backend.my_min_tmaxfrac, backend.my_max_tmaxfrac)
            

    # Obtain histogram
    pts.plot_1dhist(data, nbins,axes_labels,out_fig_filename, fig_fontsize)
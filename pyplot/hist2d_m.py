
import plottools as pts
import labels
from pylab import *
import backend

def generate_histogram_2d(ID, toplot, chainID, filenames):
    
    # Unpack filenames from input
    (priors_filename, data_filename, PLOTSDIR) = filenames
    
    # Check whether any of the given parameters are "derived",
    # since all plots containing derived parameters are appended with a character.
    isDERIVED = labels.CheckDerived(toplot, priors_filename)
    
    # How many bins for the histogram
    nbins = backend.nbins_2d

    # Font-size for the axes labels
    fig_fontsize = backend.fig_fontsize

    # Construct a filename for the output figure
    out_fig_filename = PLOTSDIR + isDERIVED + ID + '_' + toplot[0] + '_' + toplot[1] + backend.like_2d_fn

    if backend.verbose > 1:
        print 'OUTPUT TO:',out_fig_filename

    # Sort out data location & labels
    dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

    # Get the data columns
    (x,y) = pts.data_column(data_filename, dloc, fac)

    num = 0
    nt = -1
    for t in toplot:        
        if t == 'tmaxfrac':                
            nt = num
        num +=1
    
    if backend.mod_tmaxfrac:    
        (x,y) = pts.mod_data_tmax_2d((x,y), backend.my_min_tmaxfrac, backend.my_max_tmaxfrac, nt)    

    # Obtain histogram
    pts.plot_2dhist((x,y), axes_labels, nbins, fig_fontsize, out_fig_filename)
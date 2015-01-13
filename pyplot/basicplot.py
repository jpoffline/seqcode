import plottools as pts
import labels
import pylab as P
from pylab import *
data_filename = '../chains/run2/chain0001.dat'

fig_fontsize = 13
toplot_ON_x = 'm3'
toplot_ON_y = 'tmaxfrac'

out_fig_filename = toplot_ON_x +'_'+toplot_ON_y+'_corr.pdf'
dloc, lab, fac = labels.mcmc_labels((toplot_ON_x, toplot_ON_y))
data_plottable = pts.data_column(data_filename, dloc, fac)
pts.plot(data_plottable, lab, out_fig_filename, fig_fontsize)



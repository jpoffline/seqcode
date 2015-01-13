import plottools as pts
import labels
import pylab as P
from pylab import *
data_filename = '../chains/run2/chain0001.dat'
PLOTSDIR = 'plots/'

fig_fontsize = 13
toplot = []
if len(sys.argv) > 1:
    toplot.append(sys.argv[1])
else:
    toplot.append('Omk')

out_fig_filename = PLOTSDIR + toplot[0] + '_mean.pdf'

# Sort out data location & labels
dloc, axes_labels, fac = labels.mcmc_labels(toplot)

# Get the data columns
data = pts.data_1d(data_filename, dloc[0], fac[0])
ndp = len(data)

running_total = 0.0
draw_number = 1
mean = []
draw = []
for dp in xrange(0, ndp):
    running_total += data[dp]
    mean.append(running_total/draw_number)
    draw.append(draw_number)
    draw_number +=1


axes_labels.append(axes_labels[0])
axes_labels[0] = 'sample #'
pts.plot((draw,mean), axes_labels, out_fig_filename, fig_fontsize)
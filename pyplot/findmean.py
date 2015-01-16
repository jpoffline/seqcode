import plottools as pts
import labels
import pylab as P
from pylab import *


fig_fontsize = 13
toplot = []
ID = 'run'
if len(sys.argv) > 1:
    ID = sys.argv[1]
    toplot.append(sys.argv[2])
else:
    toplot.append('Omk')
    
data_DIR = '../chains/' + ID + '/'
priors_filename = data_DIR + 'priors.txt'
data_filename = data_DIR + 'chain0001.dat'
PLOTSDIR = 'plots/';    

out_fig_filename = PLOTSDIR + ID + '_' + toplot[0] + '_mean.pdf'

# Sort out data location & labels
dloc, axes_labels, fac = labels.mcmc_labels(toplot, priors_filename)

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
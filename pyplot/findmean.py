import math
import plottools as pts
import labels
import pylab as P
from pylab import *
import backend



RUN_ID = 'seq'
chainID = 'chain0001'

toplot = []
if len(sys.argv) > 1:
    RUN_ID = sys.argv[1]
    toplot.append(sys.argv[2])
else:
    toplot.append('Omk')

def ComputeMean_StandardDev(data):
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
        
    actual_mean = mean[len(mean) - 1]
    
    sd = 0.0
    for dp in xrange(0, ndp):
        sd =+ math.pow( data[dp] - actual_mean, 2.0)
    sd = math.sqrt(sd / (ndp - 1.0))    
    
    return (draw, mean, actual_mean, sd)
    
def ComputePlotMean(RUN_ID, chainID, WHATtoplot, toplot):  
      
    (priors_filename, data_filename, PLOTSDIR) = labels.GetFileNames(RUN_ID, chainID)

    # Sort out data location & labels
    dloc, axes_labels, fac = labels.mcmc_labels(WHATtoplot, priors_filename)

    # Get the data columns
    data = pts.data_1d(data_filename, dloc[0], fac[0])
    
    # Compute the mean
    (draw, mean, actual_mean, standarddev) = ComputeMean_StandardDev(data)
    
    print 'Mean of ', WHATtoplot[0], ' = ', actual_mean, 'pm', standarddev
        
    if toplot:
        DIR = backend.plotsDIR + '/' + RUN_ID + '/' + backend.meansDIR + '/'
        pts.check_dir_exists(DIR)
        out_fig_filename = DIR + RUN_ID + '_' + WHATtoplot[0] + backend.mean_fig_fn
        print 'Saving plot:', out_fig_filename
        fig_fontsize = backend.fig_fontsize
        axes_labels.append(axes_labels[0])
        axes_labels[0] = 'sample #'
        pts.plot((draw,mean), axes_labels, out_fig_filename, fig_fontsize)
    
    
    
ComputePlotMean(RUN_ID, chainID, toplot, True)    
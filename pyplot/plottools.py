from pylab import *
import pylab as P
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import os


line_BURNS = 900
tarm_loc = 8


def check_dir_exists(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        print "created output directory"

# Function to read in data (possibly with comments, marked on with '#'
# then return the required columns



def data_column(df, (p1, p2), fac):
    d = open(df,'r');
    lines = d.readlines();
    d.close();
    x = []
    y = []
    line_num = 0

    for line in lines:
    	if line.startswith('#'):
    		continue
    	else:
            if line_num > line_BURNS:
                p = line.split()
                if float(p[tarm_loc]) > 0:
                    x.append(float(p[p1]) * fac[0])
                    y.append(float(p[p2]) * fac[1])
            line_num+=1
    return (np.array(x), np.array(y))
    
def data_1d(df, p1, fac):
    d = open(df,'r');
    lines = d.readlines();
    d.close();
    x = []
    line_num = 0

    for line in lines:
    	if line.startswith('#'):
    		continue
    	else:
            if line_num > line_BURNS:
                p = line.split()
                if float(p[tarm_loc]) > 0:
                    x.append(float(p[p1]) * fac)
            line_num+=1        
    return np.array(x)    	

# Function to plot data with latex labels
def plot( (X,Y), label, output_fig_name, fig_fontsize):
    pp = PdfPages(output_fig_name)
    plt.figure()
    plt.clf()
    plt.plot(X,Y,'k',linestyle='none',marker='.')
    plt.xlabel(label[0], fontsize = fig_fontsize)
    plt.ylabel(label[1], fontsize = fig_fontsize)
    pp.savefig()
    pp.close()
    plt.close()

def plot_1dhist(data, nbins, axes_labels,out_fig_filename, fig_fontsize):
    pp = PdfPages(out_fig_filename)
    plt.figure()
    plt.clf()
    n, bins, patches = plt.hist(data, nbins, normed=1, histtype='step')
    plt.xlabel(axes_labels[0], fontsize = fig_fontsize)
    plt.ylabel('', fontsize = fig_fontsize)
    pp.savefig()
    pp.close()  
    plt.close()
    
def plot_2dhist(data_plottable, lab, nbins, fig_fontsize, output_fig_name):
    pp = PdfPages(output_fig_name)
    plt.figure()
    plt.clf()
    h,d,dd,ddd = plt.hist2d(data_plottable[0], data_plottable[1], bins = nbins)
    plt.clf()
    plt.imshow(h, aspect = "auto",origin = "lower", interpolation = "gaussian", extent = [min(data_plottable[0]),max(data_plottable[0]), min(data_plottable[1]) , max(data_plottable[1])])

    plt.xlabel(lab[0], fontsize = fig_fontsize)
    plt.ylabel(lab[1], fontsize = fig_fontsize)
    plt.colorbar()
    pp.savefig()
    pp.close() 
    plt.close()

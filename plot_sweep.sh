#!/bin/sh

####################################################
#
# Plot script: plot the sweep
#
#   J. Bloomfield (MIT) & J. Pearson (Durham)
#	(begun) March 2014
#
####################################################
#
# USAGE
#
# EXAMPLE
#
# NOTES
#
####################################################


outdir=$1
plotdir=plots
# Whats the plot going to be called?
plotname=plot_$2


# Parameter #1 column number
p1=1
# Parameter #2 column number
p2=2
# WMAP likelihood column number
WL=3
# Planck likelihood column number
PL=4
# Supernovae likelihood column number
SL=5
# Combined likelihood column number
CL=9

################################################

filename=$outdir/run$2.dat
echo $filename

gnuplot << EOF
 
set term postscript landscape color enhanced 20 solid
set output "$plotdir/$plotname.eps"

set contour base
set view map
unset key
set xlabel "w_0"
set ylabel "w_a"
unset surface

set style line 1 lc rgb "black" lw 2
set style line 2 lc rgb "magenta" lw 2
set style line 3 lc rgb "magenta" lw 2
set style line 4 lc rgb "blue" lw 2
set style line 5 lc rgb "blue" lw 2
set style line 6 lc rgb "green" lw 2
set style line 7 lc rgb "green" lw 2
set style line 8 lc rgb "red" lw 2
set style line 9 lc rgb "red" lw 2

set style increment user
set cntrparam levels discrete .63,.95
splot "$filename" u $p1:$p2:$WL w l t "WMAP",\
"" u $p1:$p2:$PL w l t "Planck",\
"" u $p1:$p2:$SL w l t "SN1a",\
"" u $p1:$p2:$CL w l t "Combined"

EOF

echo "Done plotting"
epstopdf --autorotate=All $plotdir/$plotname.eps
echo "Conversion to PDF complete"

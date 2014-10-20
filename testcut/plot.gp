set term png
set output "plot1.png"

data = "test1.dat"

set yr [2.1:3.2]


set pm3d map

set origin 0,0

set xlabel "a_cut * 1000"
set ylabel "phi_0"
set cblabel "log_10 |<R>_reg|"
splot data u ($1*1000):2:(log10(abs($3)))

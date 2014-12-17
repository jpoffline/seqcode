set terminal postscript eps enhanced color

set output "sweep_test22.eps"
d="test22.dat"
p1=1
p2=2
set xr [-3.5:-2]
p1l="log_{10}m^3"
p2l="{/Symbol W}_k h^2"

set pm3d map
set xl p1l
set yl p2l
set cbl "log_{10}|<R>|"
splot d u p1:p2:(log10(abs($4)))

reset

set output "sweep_test22_amax.eps"
d="test22.dat"
p1=1
p2=2
p1l="log_{10}m^3"
p2l="{/Symbol W}_k h^2"

set pm3d map
set xl p1l
set yl p2l
set cbl "a_{max}"
splot d u p1:p2:7


reset

set output "sweep_test22_tmax.eps"
d="test22.dat"
p1=1
p2=2
p1l="log_{10}m^3"
p2l="{/Symbol W}_k h^2"

set pm3d map
set xl p1l
set yl p2l
set cbl "t_{max}"
splot d u p1:p2:9


reset

set output "sweep_test22_w0.eps"
d="test22.dat"
p1=1
p2=2
p1l="log_{10}m^3"
p2l="{/Symbol W}_k h^2"

set pm3d map
set xl p1l
set yl p2l
set cbl "10^3(1+w_0)"
splot d u p1:p2:(1E3*(1+$8))



reset

set output "sweep_test22_w0_om.eps"
d="test22.dat"
p1=1
p2=2
p1l="log_{10}m^3"
p2l="{/Symbol W}_k h^2"

set pm3d map
set xl p1l
set yl p2l
set xr [-3.5:-2]
set cbl "(1+w_0)/m^3"
splot d u p1:p2:((1+$8)/(10**$1))



reset

set output "sweep_test22_w0_oOk.eps"
d="test22.dat"
p1=1
p2=2
p1l="log_{10}m^3"
p2l="{/Symbol W}_k h^2"

set pm3d map
set xl p1l
set yl p2l
set xr [-3.5:-2]
set yr [-0.0002:-0.001]
set cbl "(1+w_0)/{/Symbol W}_kh^2"
splot d u p1:p2:((1+$8)/$2)


#!/bin/bash
alias python=python2.7


#python findmean.py m3
#python findmean.py Omk
#python findmean.py tmaxfrac

python hist2d.py m3 phidot0
python hist2d.py m3 Omk
python hist2d.py m3 tmaxfrac
python hist2d.py Omk tmaxfrac


python hist1d.py m3  
python hist1d.py Omk
python hist1d.py tmaxfrac



import os

# Parameters to sweep over (to make plots of)
all = ('phidot0','m3','Omk','amax','tmaxfrac')

PY = 'python2.7'

for i in xrange(0, len(all)):
    print all[i]
    string_to_exe = PY + ' hist1d.py ' + all[i]
    os.system(string_to_exe)
    for j in xrange(i+1,len(all)):
        print all[i], all[j]
        string_to_exe = PY + ' hist2d.py ' + all[i] + ' ' + all[j]
        os.system(string_to_exe)    



import os, sys
import labels as lab
RUN_ID = sys.argv[1]

priors_file_name = '../chains/' + RUN_ID + '/priors.txt'

# Parameters to sweep over (to make plots of)
all = lab.read_priors(priors_file_name)
all.append('Omk')
all.append('amax')
all.append('tmaxfrac')

print all
PY = 'python2.7'

for i in xrange(0, len(all)):
    print all[i]
    string_to_exe = PY + ' hist1d.py ' + RUN_ID + ' ' + all[i]
    os.system(string_to_exe)
    for j in xrange(i+1,len(all)):
        print all[i], all[j]
        string_to_exe = PY + ' hist2d.py ' + RUN_ID + ' ' + all[i] + ' ' + all[j]
        os.system(string_to_exe)    

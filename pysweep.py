import os
import fileinput
import sys
from math import log10, pow

###################################
# BEGIN user-edit
###################################

# Min value of parameter to increment
param_min = 1.0

# Max value of parameter
param_max = 2.0

# Step size of parameter
dparam = 0.1

# Values of parameters to fix...

# phidot0 :: FIXED
param2 = 0.01
# mass :: FIXED
param3 = 1.0

###################################
# END user-edit
###################################

# now for the main guts of the sweeper code

# file name of the "prototype" parameter file (i.e., with flags)
params_proto_file = "params.proto"

# file name of the parameter file that will actually be run;
# ... this is arbitrary
params_run_file = "paramsnew.ini"

# Open up prototype parameter file
protp_file = open(params_proto_file,'r')

# read prototype into string "newparamsp"
newparamsp = protp_file.read()

#close prototype parameter file
protp_file.close()

# start param off at param_min
param = param_min

# Begin sweep
while param <= param_max:

    # increment phi0
    param1 = param
    
    # print param values to screen
    print "phi0 = " + str(param1)
    print "phidot0 = " + str(param2)
    print "mass = " + str(param3)

    # find and replace
    newparams = newparamsp.replace("param1FLAG",str(param1))
    newparams = newparams.replace("param2FLAG",str(param2))
    newparams = newparams.replace("param3FLAG",str(param3))
    
    # dump new params data to file
    p = open(params_run_file,'w')
    p.write(newparams)
    p.close()

    # run instance of deevolve with the edited .ini file
    str_run = "./main " + params_run_file
    os.system(str_run)

    # incrememnt value of the parameter
    param = param + dparam
    

print "done"

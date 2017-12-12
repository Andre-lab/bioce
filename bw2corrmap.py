"""
Reads in Bayesian fit file and calculates corrmap between experimentl and
scaled ensemble avergaed intensity.

Requires freesas to be installed: https://github.com/kif/freesas
"""
import sys
import numpy as np
from freesas.cormap import gof

fit_file = sys.argv[1]
fit_data = np.genfromtxt(fit_file)
experimental = fit_data[:,1]
simulated = fit_data[:,2]

corr = gof(experimental, simulated)
print corr

#Saving corr files
np.savetxt("exp_"+fit_file, experimental)
np.savetxt("sim_"+fit_file,simulated)

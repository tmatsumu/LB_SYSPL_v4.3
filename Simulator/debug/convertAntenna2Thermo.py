import numpy as np
import pylab as py
import sys

Tcmb = 2.725
h = 6.62606957e-34
k_b = 1.3806488e-23

nu_GHz = float(sys.argv[1])
nu = nu_GHz*1.e9

x = (h*nu)/k_b/Tcmb

print x
print np.exp(x)
print np.exp(x)-1.
print (np.exp(x)-1.)**2
print (x**2*np.exp(x))
print (np.exp(x)-1.)**2 / (x**2*np.exp(x))

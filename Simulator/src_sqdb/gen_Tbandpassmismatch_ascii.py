import numpy as np
import pylab as py
import sys

filename_out = sys.argv[1]

num = 370

mu = 1.
sigma = 0.005
alpha_d_arr = np.random.normal(mu,sigma,num)
alpha_s_arr = np.ones(num)

f = open(filename_out, "w")
for i in range(0,num):
    f.write('%d %f %f \n' % (i, alpha_d_arr[i], alpha_s_arr[i]))
print '[WRITE TEXT FILE]: ', filename_out
f.close()

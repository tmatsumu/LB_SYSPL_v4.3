import numpy as np
import pylab as py
import healpy as h
import sys

dir = '/project/projectdirs/polar/data/test_tauA/observation_20100608/scan0/'
filename = 'scan_flag.fits'
flag = h.mrdfits(dir+filename, hdu=1)
print flag

ind = np.where((flag == 1))
sys.exit()

nb = len(ind)
ind0 = ind[0]

idx_s = []
for i in range(0,nb):
    if ind[i] > ind0+1:
        idx_s.append(ind[i])
        ind0 = ind[i]

f = open(dir+'/pointing_flag.txt', "w")
f.write('%s\n' % ('scan idx, idx_start, idx_end, direction, flag'))
for i in range(0,nb_hscan):
    f.write('%d   %d   %d   %d   %d\n' % (i, idx_s[i], idx_s[i]+100, int(1), 0))
f.close()

py.plot(ind,flag[ind])
py.plot(ind_s,flag[ind_s])
py.show()

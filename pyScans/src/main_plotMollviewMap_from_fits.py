import numpy as np
import pylab as py
import healpy as h
import sys

dir_in = sys.argv[1]
filename = sys.argv[2]
max = float(sys.argv[3])
min = float(sys.argv[4])

map = h.read_map(dir_in+'/'+filename+'.fits')

max = np.max(map)
min = 0
h.mollview(map, max=max, min=min) #, nest=True)
py.savefig(dir_in+'/'+filename+'.png')
py.clf()

import numpy as np
import pylab as py
import healpy as h
import sys
import os

filename = sys.argv[1]
field = int(sys.argv[2])
max_in = float(sys.argv[3])
title = sys.argv[4]

map = h.read_map(filename, field=field)
h.mollview(map, max=max_in, min=-max_in, title=title)
py.savefig('tmp.png')
os.system('display tmp.png &')

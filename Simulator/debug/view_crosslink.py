import numpy as np
import pylab as py
import healpy as h
import sys
import os

dir = sys.argv[1]
field = int(sys.argv[2])
max_in = float(sys.argv[3])

mapA = h.read_map(dir+'/mapA_CNexcluded.fits', field=field)
mapB = h.read_map(dir+'/mapB_CNexcluded.fits', field=field)
mapH = h.read_map(dir+'/mapH.fits', field=field)

h.write_map(dir+'/mapCL_CNexcluded.fits',(mapA**2+mapB**2)/mapH**2)

h.mollview((mapA**2+mapB**2)/mapH**2, max=max_in, min=0, title='')
py.title('$<\cos{2\\alpha}>^2+<\sin{2\\alpha}>^2$')
#h.mollview((mapA**2)/mapH, max=max_in, min=0)
#h.mollview(mapA/mapH, max=max_in, min=0)
py.savefig('tmp.png')
os.system('display tmp.png &')


import numpy as np
import pylab as py
import healpy as h
import sys

dir = '/group/cmb/litebird/simdata/Maps/CMB/tmp/'
filename = 'group1_map_band1_ud512.fits'

mapT = h.read_map(dir+filename, field=0)
mapQ = h.read_map(dir+filename, field=1)
mapU = h.read_map(dir+filename, field=2)

h.mollview(mapT,title='T',max=0.0114248,min=-0.00240018)
py.savefig(dir+'mapT.png')

h.mollview(mapQ,title='Q',max=0.0001,min=-0.0001)
py.savefig(dir+'mapQ.png')

h.mollview(mapU,title='U',max=0.0001,min=-0.0001)
py.savefig(dir+'mapU.png')

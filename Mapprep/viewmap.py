import numpy as np
import pylab as py
import healpy as h
import os
import sys

filename = sys.argv[1]
TQU = sys.argv[2]

if TQU == 'T':
    map = h.read_map(filename+'.fits')
    h.mollview(map)
    py.savefig(filename+'.png')
    os.system('display '+filename+'.png')


if TQU == 'TQU':
    map = h.read_map(filename+'.fits',field=0)
    h.mollview(map)
    py.savefig(filename+'_T.png')
    os.system('display '+filename+'_T.png')

    map = h.read_map(filename+'.fits',field=1)
    h.mollview(map)
    py.savefig(filename+'_Q.png')
    os.system('display '+filename+'_Q.png')

    map = h.read_map(filename+'.fits',field=2)
    h.mollview(map)
    py.savefig(filename+'_U.png')
    os.system('display '+filename+'_U.png')

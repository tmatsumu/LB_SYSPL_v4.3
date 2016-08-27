import numpy as np
import pylab as py
import healpy as h
import sys

pi = np.pi
radeg = (180./pi)

dir = sys.argv[1]
filename = (sys.argv[2])
nsideout = int(sys.argv[3])
beam_size = int(sys.argv[4])
dirout = sys.argv[5]

T = h.read_map(dir+filename+'.fits',field=0)
T = h.ud_grade(T, nsideout)
Q = h.read_map(dir+filename+'.fits',field=1)
Q = h.ud_grade(Q, nsideout)
U = h.read_map(dir+filename+'.fits',field=2)
U = h.ud_grade(U, nsideout)

if beam_size != 0:
    print '[gen_mapprep.py] smoothing to ', beam_size, 'arcmin'
    T = h.smoothing(T,fwhm=beam_size/60./radeg)
    Q = h.smoothing(Q,fwhm=beam_size/60./radeg)
    U = h.smoothing(U,fwhm=beam_size/60./radeg)

TQU = [T,Q,U]
h.write_map(dirout+filename+'_ud'+str(nsideout)+'_b'+str(beam_size)+'arcmin.fits',TQU)
h.write_map(dirout+filename+'_T_ud'+str(nsideout)+'_b'+str(beam_size)+'arcmin.fits',T)


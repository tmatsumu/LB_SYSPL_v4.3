import numpy as np
import pylab as py
import healpy as h
import sys
import os

pi = np.pi
radeg = (180./pi)

def gen_dipole(nside,unit):
    # from WMAP 1st yaer Bennett et al
    l = 263.85 # +/-0.1 degs
    b = 48.25 # +/-0.04 degs
    amp_mK = 3.346 # +/-0.017 mK
    if unit == 'uK': amp_unit = amp_mK * 1.e3; print unit
    if unit == 'mK': amp_unit = amp_mK; print unit
    if unit == 'K': amp_unit = amp_mK * 1.e-3; print unit

    npix = h.nside2npix(nside)
    ipix = range(npix)
    theta, phi = h.pix2ang(nside,ipix)
    dipole = np.cos(theta)
    
    xyz = h.ang2vec(theta,phi)
    theta0 = pi/2.-b/radeg;  phi0 = l/radeg
    x = np.cos(phi0)*np.cos(theta0)*xyz[:,0] + np.sin(phi0)*np.cos(theta0)*xyz[:,1] - np.sin(theta0)*xyz[:,2]
    y = - np.sin(phi0)*xyz[:,0]+np.cos(phi0)*xyz[:,1]
    z = np.sin(theta0)*np.cos(phi0)*xyz[:,0] + np.sin(phi0)*np.sin(theta0)*xyz[:,1] + np.cos(theta0)*xyz[:,2]

    ipix = h.vec2pix(nside,x,y,z)

    dipole_out = amp_unit*dipole[ipix]
    return dipole_out

#++++++++++++++++++++++++++++++++++++++++

fname_map = sys.argv[1]
print fname_map

mapT = h.read_map(fname_map+'.fits',field=0)
mapQ = h.read_map(fname_map+'.fits',field=1)
mapU = h.read_map(fname_map+'.fits',field=2)

npix = len(mapT)
nside = h.npix2nside(npix)

dipole = gen_dipole(nside,'K')
#h.mollview(dipole)
#py.savefig('dipole.png')
#os.system('display dipole.png&')
#sys.exit()

mapTQU = [mapT*1.e-6+dipole,mapQ*1.e-6,mapU*1.e-6]
#mapTQU = [mapT,mapQ,mapU]
h.write_map(fname_map+'_dipole.fits',mapTQU)
h.mollview(mapTQU[0])

py.savefig(fname_map+'_dipole.png')
#os.system('display '+fname_map+'_dipole.png')
print ''
